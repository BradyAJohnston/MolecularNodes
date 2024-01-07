from pathlib import Path

from ...blender import obj, coll
from .molecule import Molecule
from .bcif import BCIF
from .pdbx import PDBX

class Ensemble(Molecule):
    def __init__(self, file_path):
        self.file_path = file_path
        self.molecule = self.read_file(self.file_path)
        self.file = self.molecule.file
    
    def read_file(self):
        ext = Path(self.file_path).suffix
        if ext in (".bcif", ".bin"):
            data = BCIF(self.file_path)
        else:
            data = PDBX(self.file_path, extra_fields = ['label_entity_id'])
        
        return data
    
    def create_data_object(self, name = 'NewEnsemble'):
        obj_data = obj.create_data_object(self.assemblies(as_array=True), name = name, collection = coll.mn())
        return obj_data