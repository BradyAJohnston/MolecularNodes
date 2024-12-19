from pathlib import Path

import numpy as np
import bpy

from .ensemble import Ensemble
from .bcif import BCIF
from .cif import OldCIF
from ..molecule import molecule
from ... import blender as bl
from databpy import store_named_attribute, AttributeTypes
from ... import color


class CellPack(Ensemble):
    def __init__(self, file_path):
        super().__init__(file_path)
        self.file_type = self._file_type()
        self.chain_ids: np.ndarray
        self.transformations: np.ndarray
        self.array = self._read(self.file_path)

    def create_object(
        self,
        name="CellPack",
        node_setup: bool = True,
        world_scale: float = 0.01,
        fraction: float = 1.0,
        simplify=False,
    ):
        self.object = self._create_data_object(name=f"{name}")
        self._create_object_instances(name=name, node_setup=node_setup)
        self._setup_node_tree(fraction=fraction)
        return self.object

    def _file_type(self):
        return Path(self.file_path).suffix.strip(".")

    def _read(self, file_path):
        "Read a Cellpack File"
        suffix = Path(file_path).suffix

        if suffix in (".bin", ".bcif"):
            data = BCIF(file_path)
        elif suffix == ".cif":
            data = OldCIF(file_path)
        else:
            raise ValueError(f"Invalid file format: '{suffix}")

        self.chain_ids = data.chain_ids
        self.transformations = data.assemblies(as_array=True)
        return data.array

    def _create_object_instances(
        self, name: str = "CellPack", node_setup: bool = True
    ) -> bpy.types.Collection:
        collection = bl.coll.cellpack(name)

        if self.file_type == "cif":
            array = self.array[0]
        else:
            array = self.array
        for i, chain in enumerate(np.unique(array.chain_id)):
            chain_atoms = array[array.chain_id == chain]
            obj_name = f"{str(i).rjust(4, '0')}_{chain}"

            obj, coll_none = molecule._create_object(
                array=chain_atoms,
                name=obj_name,
                collection=collection,
            )

            colors = np.tile(color.random_rgb(i), (len(chain_atoms), 1))
            store_named_attribute(
                obj=obj,
                data=colors,
                name="Color",
                atype=AttributeTypes.FLOAT_COLOR,
            )

            if node_setup:
                bl.nodes.create_starting_node_tree(
                    obj, name=f"MN_pack_instance_{name}", color=None
                )

        self.data_collection = collection
        self.instance_collection = collection

        return collection

    def _create_data_object(self, name="DataObject"):
        data_object = bl.mesh.create_data_object(
            self.transformations, name=name, collection=bl.coll.mn()
        )

        data_object["chain_ids"] = self.chain_ids

        return data_object

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
