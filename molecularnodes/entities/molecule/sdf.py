from .molecule import Molecule
from biotite.structure.io.mol import MOLFile


class SDF(Molecule):
    def __init__(self, file_path):
        super().__init__(file_path)
        self.file = self.read(self.file_path)
        self.array = self._get_structure()
        self.n_atoms = self.array.array_length()

    def read(self, file_path):
        return MOLFile.read(file_path)

    def _get_structure(self):
        return self.file.get_structure()

    def _assemblies(self):
        # TODO maybe look into symmetry operations for small mols
        return None
