from abc import ABCMeta
from pathlib import Path
from io import BytesIO
from biotite.structure import AtomArray, AtomArrayStack
from biotite import structure as struc
import numpy as np
from ...assets import data
from ... import color


class ReaderBase(metaclass=ABCMeta):
    """
    Abstract base class for reading molecular data from a file.
    """

    def __init__(self, file_path: str | Path | BytesIO):
        self.file_path = file_path
        self.file = self.read(file_path)
        self.array = self.get_structure()
        self.set_standard_annotations()

    def read(self, file_path):
        raise NotImplementedError("Subclasses must implement this method.")

    def set_standard_annotations(self):
        annotations = {
            "mass": _compute_mass,
            "atomic_number": _compute_atomic_number,
            "res_name_int": _compute_res_name_int,
            "chain_id_int": _compute_chain_id_int,
            "vdw_radii": _compute_vdw_radii,
            "atom_name_int": _compute_atom_name_int,
            "charge": _compute_charge,
            "lipophobicity": _compute_lipophobicity,
            "Color": _compute_color,
            "is_alpha_carbon": _compute_is_alpha_carbon,
            "is_solvent": _compute_is_solvent,
            "is_backbone": _compute_is_backbone,
            "is_nucleic": _compute_is_nucleic,
            "is_peptide": _compute_is_peptide,
            "is_side_chain": _compute_is_side_chain,
            "is_carb": _compute_is_carb,
        }
        for key, func in annotations.items():
            try:
                self.array.set_annotation(key, func(self.array))
            except Exception:
                pass

    def get_structure(self) -> AtomArrayStack | AtomArray:
        raise NotImplementedError("Subclasses must implement this method.")


def _compute_mass(array):
    return np.array(
        [
            data.elements.get(x, {}).get("standard_mass", 0.0)
            for x in np.char.title(array.element)
        ]
    )


def _compute_atomic_number(array):
    return np.array(
        [
            data.elements.get(x, {}).get("atomic_number", 0)
            for x in np.char.title(array.element)
        ]
    )


def _compute_res_name_int(array):
    return np.array(
        [
            data.residues.get(name, {}).get("res_name_num", -1)
            for name in array.res_name
        ],
        dtype=int,
    )


def _compute_chain_id_int(array):
    if isinstance(array.chain_id[0], int):
        return array.chain_id
    else:
        return np.unique(array.chain_id, return_inverse=True)[1]


def _compute_vdw_radii(array):
    return np.array(
        [
            data.elements.get(x, {}).get("vdw_radii", 100.0) / 100
            for x in np.char.title(array.element)
        ]
    )


def _compute_atom_name_int(array):
    return np.array([data.atom_names.get(x, -1) for x in array.atom_name])


def _compute_charge(array):
    return np.array(
        [
            data.atom_charge.get(res, {}).get(atom, 0)
            for res, atom in zip(array.res_name, array.atom_name)
        ]
    )


def _compute_lipophobicity(array):
    return np.array(
        [
            data.lipophobicity.get(res, {}).get(atom, 0)
            for res, atom in zip(array.res_name, array.atom_name)
        ]
    )


def _compute_color(array, color_plddt: bool = False):
    if color_plddt:
        return color.plddt(array.b_factor)
    else:
        return color.color_chains(array.atomic_number, array.chain_id)


def _compute_is_alpha_carbon(array):
    return array.atom_name == "CA"


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


def _compute_is_solvent(array):
    return struc.filter_solvent(array)


def _compute_is_nucleic(array):
    return struc.filter_nucleotides(array)


def _compute_is_peptide(array):
    return np.logical_or(
        struc.filter_amino_acids(array),
        struc.filter_canonical_amino_acids(array),
    )


def _compute_is_side_chain(array):
    not_backbone = np.logical_not(array.is_backbone)

    return np.logical_and(
        not_backbone, np.logical_or(array.is_nucleic, array.is_peptide)
    )


def _compute_is_carb(array):
    return struc.filter_carbohydrates(array)
