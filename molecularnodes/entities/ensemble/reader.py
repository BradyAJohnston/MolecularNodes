from io import BytesIO
from .bcif import BCIF
from ..molecule.pdbx import PDBX
from biotite.structure.io import pdbx
from biotite import structure as struc
from biotite import InvalidFileError
from pathlib import Path
import numpy as np


# we have to override the reading of lines to strip out the preceeding whitespace
# for the petworld files
class PetworldCIFFileReader(pdbx.CIFFile):
    @classmethod
    def read(cls, file):
        """
        Read a CIF file.

        Parameters
        ----------
        file : file-like object or str
            The file to be read.
            Alternatively a file path can be supplied.

        Returns
        -------
        file_object : CIFFile
            The parsed file.
        """
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
        self._extra_annotations["asym_id"] = self._set_asym_id
        self.file_path = file_path
        self.file: pdbx.BinaryCIFFile | pdbx.CIFFile = self._read()
        self.n_molecules: int = pdbx.get_model_count(self.file)
        self.molecules: list[stuc.AtomArray] = []

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

    def _set_asym_id(self, array, file) -> np.ndarray:
        return array.chain_id

    def get_structure(
        self,
        extra_fields=["b_factor", "occupancy", "atom_id"],
        bonds=True,
        model: int | None = None,
    ):
        array: struc.AtomArray = pdbx.get_structure(
            self.file, model=model, extra_fields=extra_fields
        )  # type: ignore
        array = self.set_extra_annotations(array, self.file)

        if not array.bonds and bonds:
            array.bonds = struc.bonds.connect_via_residue_names(
                array, inter_residue=True
            )

        print(f"{array[:10]=}")
        return array

    def get_molecules(
        self, extra_fields=["b_factor", "occupancy", "atom_id"], bonds=True
    ):
        try:
            array = self.get_structure(extra_fields, bonds)
            self.molecules = list(
                [array[array.chain_id == c] for c in np.unique(array.chain_id)]
            )
        except InvalidFileError as e:
            for i in range(self.n_molecules):
                self.molecules.append(self.get_structure(model=int(i + 1)))
