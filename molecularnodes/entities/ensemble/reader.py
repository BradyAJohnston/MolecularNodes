from .bcif import BCIF
from ..molecule.pdbx import PDBX
from biotite.structure.io import pdbx
from pathlib import Path


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
        self.file_path = file_path
        self.file = self._read()
        self.array = self.get_structure()

    def _read(self):
        suffix = Path(self.file_path).suffix
        if suffix == ".bcif":
            return pdbx.BinaryCIFFile.read(self.file_path)
        elif suffix == ".cif":
            return PetworldCIFFileReader.read(self.file_path)
        else:
            raise ValueError(f"Invalid file format: '{suffix}")
