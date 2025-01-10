from io import BytesIO
from pathlib import Path

import numpy as np
from biotite import InvalidFileError
from biotite import structure as struc
from biotite.structure.io import pdbx

from ..molecule.pdbx import PDBX


# For reading cellpack files, we override the CIFFile from biotite. The only change we
# implement is that we add the `line.strip()` when initially getting the lines from the
# text in the `deserialize` method. This fixes the reading of the PETWORLD files, and we
# don't have to write-out a modifier version of the file before reading it back in
class PetworldCIFFileReader(pdbx.CIFFile):
    @classmethod
    def read(cls, file):
        # File name
        if pdbx.cif.is_open_compatible(file):
            with open(file, "r") as f:
                text = f.read()
        # File object
        else:
            if not pdbx.cif.is_text(file):
                raise TypeError("A file opened in 'text' mode is required")
            text = file.read()
        return cls.deserialize(text)

    @staticmethod
    def deserialize(text):
        lines = [line.strip() for line in text.splitlines()]
        block_starts = []
        block_names = []
        for i, line in enumerate(lines):
            if not pdbx.cif._is_empty(line):
                data_block_name = pdbx.cif._parse_data_block_name(line)
                if data_block_name is not None:
                    block_starts.append(i)
                    block_names.append(data_block_name)
        return pdbx.CIFFile(
            pdbx.cif._create_element_dict(lines, block_names, block_starts)
        )


class CellPackReader(PDBX):
    def __init__(self, file_path):
        super().__init__(file_path)
        self._extra_annotations["asym_id"] = self._get_asym_id
        self._extra_annotations["pdb_model_num"] = self._get_pdb_model_num
        self.file_path = file_path
        self.file: pdbx.BinaryCIFFile | pdbx.CIFFile = self._read()
        self.n_molecules: int = pdbx.get_model_count(self.file)
        self.molecules: dict[str, stuc.AtomArray] = {}

    @property
    def mol_ids(self) -> np.ndarray:
        return np.unique(list(self.molecules.keys()))

    def _parse_filepath(self, file_path: Path | str | BytesIO) -> None:
        pass

    def _read(self):
        suffix = Path(self.file_path).suffix
        if suffix == ".bcif":
            return pdbx.BinaryCIFFile.read(self.file_path)
        elif suffix == ".cif":
            return PetworldCIFFileReader.read(self.file_path)
        else:
            raise ValueError(f"Invalid file format: '{suffix}")

    def _get_asym_id(self, array, file) -> np.ndarray:
        return array.chain_id

    def _get_pdb_model_num(self, array, file) -> np.ndarray:
        return (
            list(self.file.values())[0]["atom_site"]["pdbx_PDB_model_num"]
            .as_array()
            .astype(int)
        )

    def get_structure(
        self,
        extra_fields=["b_factor", "occupancy", "atom_id"],
        bonds=True,
        model: int | None = None,
    ) -> struc.AtomArray:
        array: struc.AtomArray = pdbx.get_structure(
            self.file, model=model, extra_fields=extra_fields
        )  # type: ignore
        array = self.set_extra_annotations(array, self.file)

        if not array.bonds and bonds:
            array.bonds = struc.bonds.connect_via_residue_names(
                array, inter_residue=True
            )
        return array

    @property
    def blocks(self):
        return list(self.file.values())[0]

    def get_molecules(
        self, extra_fields=["b_factor", "occupancy", "atom_id"], bonds=True
    ):
        self._is_petworld = False

        if "PDB_model_num" in self.blocks["pdbx_struct_assembly_gen"]:
            self._is_petworld = True

        try:
            array = self.get_structure(extra_fields, bonds)
            if isinstance(array, struc.AtomArrayStack):
                array = array[0]

            # if self._is_petworld:
            #     array.set_annotation(
            #         "chain_id", np.char.rjust(array.pdb_model_num, 4, "0")
            #     )

            self.molecules = {
                c: array[array.chain_id == c] for c in np.unique(array.chain_id)
            }
        except InvalidFileError:
            self._is_petworld = True
            for i in range(self.n_molecules):
                array = self.get_structure(extra_fields, bonds, model=i + 1)
                array.set_annotation("pdbx_PDB_model_num", np.repeat(i + 1, len(array)))
                array.chain_id = array.pdbx_PDB_model_num
                chain_name = "{}_{}".format(
                    str(i).rjust(4, "0"), str(array.chain_id[0])
                )
                self.molecules[chain_name] = array
