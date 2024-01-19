from .molecule import Molecule


class SDF(Molecule):
    def __init__(self, file_path):
        self.file_path = file_path
        self.file = self.read(self.file_path)
        self.array = self._get_structure()
        self.n_models = self.array.shape[0]
        self.n_atoms = self.array.array_length()

    def read(self, file_path):
        from biotite.structure.io.mol import MOLFile

        return MOLFile.read(file_path)

    def _get_structure(self):
        return self.file.get_structure()

    def _assemblies(self):
        # TODO maybe look into symmetry operations for small mols
        return None
