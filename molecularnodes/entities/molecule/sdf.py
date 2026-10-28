from biotite.structure import AtomArray
from biotite.structure.io.mol import MOLFile
from .reader import ReaderBase


class SDFReader(ReaderBase[MOLFile]):
    _extra_annotations = {}

    def read(self, file_path) -> MOLFile:
        return MOLFile.read(file_path)

    def get_structure(self, model: int | None = None) -> AtomArray:
        return self.file.get_structure()

    def _assemblies(self):
        # TODO maybe look into symmetry operations for small mols
        return None
