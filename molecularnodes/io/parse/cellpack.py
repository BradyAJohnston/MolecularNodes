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
        self.data = self._read(self.file_path)
        self.transformations = self._parse_transformations(self.data)
        self.chain_ids = np.unique(self.data.structure.chain_id)
    
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

    def _parse_transformations(self, data: Molecule) -> np.ndarray:
        transformations = data.assemblies
        if isinstance(transformations, dict):
            transformations = utils.array_quaternions_from_dict(transformations)
        
        return transformations
    
    def _create_object_instances(
        self,
        # array, 
        name: str = 'CellPack', 
        node_setup: bool = True
        ) -> bpy.types.Collection:
        import biotite.structure as struc
        
        if isinstance(array, struc.AtomArrayStack):
            array = array[0]

        collection = coll.cellpack(name)
        
        for i, chain in enumerate(self.chain_ids):
            # model = molecule._create_model(
            #     array = self.structure[self.structure.chain_id == chain],
            #     name=name,
            #     collection=collection
            # )
            atoms = array[array.chain_id == chain]
            model, collection = load.create_model(
                array=atoms,
                name=f"{str(i).rjust(4, '0')}_{chain}",
                collection=collection
            )
        
            # color each created model differently
            # currently this is done randomly, but should be able to support palettes without
            # too much trouble
            colors = np.tile(color.random_rgb(i), (len(atoms), 1))
            obj.add_attribute(model, name="Color", data=colors, type="FLOAT_COLOR", overwrite=True)
            
            if node_setup:
                nodes.create_starting_node_tree(model, name = f"MN_pack_instance_{name}", set_color=False)
        
        return collection