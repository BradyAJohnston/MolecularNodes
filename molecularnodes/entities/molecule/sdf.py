from biotite.structure.io.mol import MOLFile
from .reader import ReaderBase


class SDFReader(ReaderBase):
    _extra_annotations = {}

    def read(self, file_path):
        return MOLFile.read(file_path)

    def get_structure(self):
        return self.file.get_structure()

    def _assemblies(self):
        # TODO maybe look into symmetry operations for small mols
        return None
