from pathlib import Path
import numpy as np
import bpy

from .molecule import Molecule
from .ensemble import Ensemble
from .bcif import BCIF
from .pdbx import PDBX
from ..parse import molecule
from ...blender import coll, obj, nodes
from ... import utils, color

class CellPack(Ensemble):
    def __init__(self, file_path):
        self.file_path = file_path
        self.file_type = self._file_type()
        self.data = self._read(self.file_path)
        self.structure = self.data.structure
        self.transformations = self._parse_transformations()
        self.chain_ids = np.unique(self.data.structure.chain_id)
        
    
    def _file_type(self):
        return Path(self.file_path).suffix.strip(".")
    
    def _read(self, file_path):
        "Read a Cellpack File"
        suffix = Path(file_path).suffix
        
        if suffix in (".bin", ".bcif"):
            data = BCIF(file_path)
        elif suffix == ".cif":
            data = PDBX(file_path)
        else:
            raise ValueError(f"Invalid file format: '{suffix}")
        
        return data

    def _parse_transformations(self) -> np.ndarray:
        # TODO this is ugly as cif assemblies is a funciton and bcif returns the array
        # I need to equalise the two for consistency
        if self.file_type == "cif":
            return utils.array_quaternions_from_dict(self.data.assemblies())
        else:
            return self.data.assemblies
    
    def _create_object_instances(
        self,
        name: str = 'CellPack', 
        node_setup: bool = True
        ) -> bpy.types.Collection:

        collection = coll.cellpack(name)
        
        if self.file_type == "cif":
            array = self.structure[0]
        else:
            array = self.structure
        for i, chain in enumerate(np.unique(array.chain_id)):
            chain_atoms = array[array.chain_id == chain]
            model, coll_none = molecule._create_model(
                array = chain_atoms,
                name=f"{str(i).rjust(4, '0')}_{chain}",
                collection=collection
            )
        
            # color each created model differently
            # currently this is done randomly, but should be able to support palettes without
            # too much trouble
            colors = np.tile(color.random_rgb(i), (len(chain_atoms), 1))
            obj.add_attribute(model, name="Color", data=colors, type="FLOAT_COLOR", overwrite=True)
            
            if node_setup:
                nodes.create_starting_node_tree(model, name = f"MN_pack_instance_{name}", set_color=False)
        
        return collection