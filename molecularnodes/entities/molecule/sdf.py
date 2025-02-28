from .reader import ReaderBase
from biotite.structure.io.mol import MOLFile


class SDFReader(ReaderBase):
    def __init__(self, file_path):
        super().__init__(file_path)

    def read(self, file_path):
        return MOLFile.read(file_path)

    def get_structure(self):
        return self.file.get_structure()

    def _assemblies(self):
        # TODO maybe look into symmetry operations for small mols
        return None
