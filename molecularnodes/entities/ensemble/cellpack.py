from pathlib import Path

import numpy as np
import bpy

from .base import Ensemble
from .bcif import BCIF
from .cif import OldCIF
from ..molecule import molecule
from ... import blender as bl
from databpy import BlenderObject
from databpy import store_named_attribute, AttributeTypes
from ... import color
from .reader import CellPackReader


class CellPack(Ensemble):
    def __init__(self, file_path):
        super().__init__(file_path)
        self.file_type = self._file_type()
        self.file = CellPackReader(file_path)
        self.file.get_molecules()
        self.transformations = self.file.assemblies(as_array=True)

    @property
    def molecules(self):
        return self.file.molecules

    def create_object(
        self,
        name="CellPack",
        node_setup: bool = True,
        world_scale: float = 0.01,
        fraction: float = 1.0,
        simplify=False,
    ):
        self.object = self._create_data_object(name=f"{name}")
        self._create_object_instances(name=self.object.name, node_setup=node_setup)
        self._setup_node_tree(fraction=fraction)
        return self.object

    def _file_type(self):
        return Path(self.file_path).suffix.strip(".")

    def _create_object_instances(
        self, name: str = "CellPack", node_setup: bool = True
    ) -> bpy.types.Collection:
        collection = bl.coll.cellpack(name)

        for i, mol_id in enumerate(self.file.mol_ids):
            array = self.molecules[mol_id]
            chain_name = array.asym_id[0]

            obj, coll_none = molecule._create_object(
                array=array,
                name=mol_id,
                collection=collection,
            )

            colors = np.tile(color.random_rgb(i), (len(array), 1))
            store_named_attribute(
                obj=obj,
                data=colors,
                name="Color",
                atype=AttributeTypes.FLOAT_COLOR,
            )

            if node_setup:
                bl.nodes.create_starting_node_tree(
                    obj,
                    name=f"MN_pack_instance_{name}",
                    color=None,
                    material="MN Ambient Occlusion",
                )

        self.data_collection = collection
        self.instance_collection = collection

        return collection

    def _create_data_object(self, name="DataObject"):
        bob = BlenderObject(
            bl.mesh.create_data_object(
                self.transformations, name=name, collection=bl.coll.mn()
            )
        )
        bob.object["chain_ids"] = self.file.mol_ids

        # if we are dealing with petworld data, overwrite the chain_id for the data object
        if self.file._is_petworld:
            bob.store_named_attribute(bob.named_attribute("pdb_model_num"), "chain_id")

        return bob.object

    def _setup_node_tree(self, name="CellPack", fraction=1.0, as_points=False):
        mod = bl.nodes.get_mod(self.object)

        group = bl.nodes.new_tree(name=f"MN_ensemble_{name}", fallback=False)
        mod.node_group = group

        node_pack = bl.nodes.add_custom(group, "Ensemble Instance", location=[-100, 0])
        node_pack.inputs["Instances"].default_value = self.data_collection
        node_pack.inputs["Fraction"].default_value = fraction
        node_pack.inputs["As Points"].default_value = as_points

        link = group.links.new
        link(bl.nodes.get_input(group).outputs[0], node_pack.inputs[0])
        link(node_pack.outputs[0], bl.nodes.get_output(group).inputs[0])
