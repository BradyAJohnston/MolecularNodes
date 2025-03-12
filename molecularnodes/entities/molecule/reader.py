from abc import ABCMeta
from pathlib import Path
from io import BytesIO
from biotite.structure import AtomArray, AtomArrayStack
from biotite import structure as struc
from biotite.file import File, InvalidFileError
import numpy as np
import json
from ...assets import data
from ... import color


class ReaderBase(metaclass=ABCMeta):
    """
    Abstract base class for reading molecular data from a file.
    """

    def __init__(self, file_path: str | Path | BytesIO):
        self._extra_annotations: dict
        self.file_path = file_path
        self.file = self.read(file_path)
        self.array = self.get_structure()
        self.array = self.set_extra_annotations(
            self._extra_annotations, self.array, self.file
        )
        self.array = self.set_standard_annotations(self.array)

    @property
    def n_models(self) -> int:
        return self.array.stack_depth

    def read(self, file_path: str | Path | BytesIO) -> File:
        raise NotImplementedError("Subclasses must implement this method.")

    @classmethod
    def set_standard_annotations(cls, array):
        annotations = {
            "mass": cls._compute_mass,
            "atomic_number": cls._compute_atomic_number,
            "res_name_int": cls._compute_res_name_int,
            "chain_id_int": cls._compute_chain_id_int,
            "vdw_radii": cls._compute_vdw_radii,
            "atom_name_int": cls._compute_atom_name_int,
            "charge": cls._compute_charge,
            "lipophobicity": cls._compute_lipophobicity,
            "Color": cls._compute_color,
            "is_alpha_carbon": cls._compute_is_alpha_carbon,
            "is_solvent": cls._compute_is_solvent,
            "is_backbone": cls._compute_is_backbone,
            "is_nucleic": cls._compute_is_nucleic,
            "is_peptide": cls._compute_is_peptide,
            "is_hetero": cls._compute_is_hetero,
            "is_side_chain": cls._compute_is_side_chain,
            "is_carb": cls._compute_is_carb,
        }
        for key, func in annotations.items():
            try:
                array.set_annotation(key, func(array))
            except Exception:
                pass

        return array

    def get_structure(self, model: int | None = None) -> AtomArrayStack | AtomArray:
        raise NotImplementedError("Subclasses must implement this method.")

    def entity_ids(self) -> list[str] | None:
        return None

    def chain_ids(self) -> list[str] | None:
        return np.unique(self.array.chain_id).tolist()

    def assemblies(self, as_json_string: bool = False) -> dict | str:
        try:
            if as_json_string:
                return json.dumps(self._assemblies())
            return self._assemblies()
        except InvalidFileError:
            return ""

    def _assemblies(self):
        return {}

    @staticmethod
    def set_extra_annotations(
        annotations: dict, array, file
    ) -> AtomArray | AtomArrayStack:
        for name, func in annotations.items():
            try:
                # for the getting of some custom attributes, we get it for the full atom_site
                # but we need to assign a subset of the array that is that model in the full
                # atom_site, so we subset using the atom_id
                try:
                    array.set_annotation(name, func(array, file))
                except IndexError:
                    array.set_annotation(
                        name,
                        func(array, file)[array.atom_id - 1],  # type: ignore
                    )
            except KeyError:
                pass
                # if True:
                #     print(f"Unable to add {name} as an attribute, error: {e}")

        return array

    @staticmethod
    def _compute_mass(array):
        return np.array(
            [
                data.elements.get(x, {}).get("standard_mass", 0.0)
                for x in np.char.title(array.element)
            ],
            dtype=float,
        )

    @staticmethod
    def _compute_atomic_number(array):
        return np.array(
            [
                data.elements.get(x, {}).get("atomic_number", 0)
                for x in np.char.title(array.element)
            ],
            dtype=int,
        )

    @staticmethod
    def _compute_res_name_int(array):
        return np.array(
            [
                data.residues.get(name, {}).get("res_name_num", -1)
                for name in array.res_name
            ],
            dtype=int,
        )

    @staticmethod
    def _compute_chain_id_int(array):
        if isinstance(array.chain_id[0], int):
            return array.chain_id.astype(int)
        else:
            return np.unique(array.chain_id, return_inverse=True)[1].astype(int)

    @staticmethod
    def _compute_vdw_radii(array):
        return np.array(
            [
                data.elements.get(x, {}).get("vdw_radii", 100.0) / 100
                for x in np.char.title(array.element)
            ],
            dtype=float,
        )

    @staticmethod
    def _compute_atom_name_int(array):
        return np.array(
            [data.atom_names.get(x, -1) for x in array.atom_name], dtype=int
        )

    @staticmethod
    def _compute_charge(array):
        return np.array(
            [
                data.atom_charge.get(res, {}).get(atom, 0)
                for res, atom in zip(array.res_name, array.atom_name)
            ],
            dtype=float,
        )

    @staticmethod
    def _compute_lipophobicity(array):
        return np.array(
            [
                data.lipophobicity.get(res, {}).get(atom, 0)
                for res, atom in zip(array.res_name, array.atom_name)
            ],
            dtype=float,
        )

    @staticmethod
    def _compute_color(array, color_plddt: bool = False):
        if color_plddt:
            return color.plddt(array.b_factor)
        else:
            return color.color_chains(array.atomic_number, array.chain_id)

    @staticmethod
    def _compute_is_alpha_carbon(array):
        return array.atom_name == "CA"

    @staticmethod
    def _compute_is_hetero(array):
        return array.hetero

    @staticmethod
    def _compute_is_backbone(array):
        """
        Get the atoms that appear in peptide backbone or nucleic acid phosphate backbones.
        Filter differs from the Biotite's `struc.filter_peptide_backbone()` in that this
        includes the peptide backbone oxygen atom, which biotite excludes. Additionally
        this selection also includes all of the atoms from the ribose in nucleic acids,
        and the other phosphate oxygens.
        """
        backbone_atom_names = [
            # Peptide backbone atoms
            "N",
            "C",
            "CA",
            "H",
            "HA",
            "O",
            # Continuous nucleic backbone atoms
            "P",
            "O5'",
            "C5'",
            "C4'",
            "C3'",
            "O3'",
            # Alternative names for phosphate O's
            "O1P",
            "OP1",
            "O2P",
            "OP2",
            # Remaining ribose atoms
            "O4'",
            "C1'",
            "C2'",
            "O2'",
        ]

        is_backbone_atom = np.isin(array.atom_name, backbone_atom_names)
        is_not_solvent = np.logical_not(struc.filter_solvent(array))

        return np.logical_and(is_backbone_atom, is_not_solvent)

    @staticmethod
    def _compute_is_solvent(array):
        return struc.filter_solvent(array)

    @staticmethod
    def _compute_is_nucleic(array):
        return struc.filter_nucleotides(array)

    @staticmethod
    def _compute_is_peptide(array):
        return np.logical_or(
            struc.filter_amino_acids(array),
            struc.filter_canonical_amino_acids(array),
        )

    @staticmethod
    def _compute_is_side_chain(array):
        not_backbone = np.logical_not(array.is_backbone)

        return np.logical_and(
            not_backbone, np.logical_or(array.is_nucleic, array.is_peptide)
        )

    @staticmethod
    def _compute_is_carb(array):
        return struc.filter_carbohydrates(array)
