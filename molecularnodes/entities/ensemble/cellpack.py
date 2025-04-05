import json
from pathlib import Path
import bpy
import numpy as np
from biotite.structure import AtomArray
from databpy import AttributeTypes, BlenderObject, store_named_attribute
from ... import blender as bl
from ... import color
from ...nodes import nodes
from ..utilities import create_object
from .base import Ensemble
from .reader import CellPackReader


class CellPack(Ensemble):
    def __init__(self, file_path):
        super().__init__(file_path)
        self.file_type = self._file_type()
        self.file = CellPackReader(file_path)
        self.file.get_molecules()
        self.transformations = self.file.get_assemblies()
        self.color_entity = {}
        self._color_palette_path = Path(file_path).parent / "color_palette.json"
        self.object = self._create_data_object(name=f"{Path(file_path).name}")
        self._create_object_instances(name=self.object.name, node_setup=False)
        self._setup_node_tree(fraction=0.1)
        # self._setup_colors()

    def _setup_colors(self):
        if self._color_palette_path.exists():
            self.color_palette = json.load(open(self._color_palette_path))

        for entity in np.unique(self.array.entity_id):
            ename = self.data.entities[entity]
            if ename in self.color_palette:
                rgb = [self.color_palette[ename][c] / 255.0 for c in "xyz"]
                self.color_entity[entity] = np.array([*rgb, 1.0])
            else:
                self.color_entity[entity] = color.random_rgb(int(entity))

            self.entity_chains[entity] = (
                np.unique(self.array.asym_id[self.array.entity_id == entity]) @ property
            )

    @property
    def molecules(self):
        return self.file.molecules

    def _file_type(self):
        return Path(self.file_path).suffix.strip(".")

    def _assign_colors(self, obj: bpy.types.Object, array: AtomArray):
        # random color per chain
        # could also do by entity, + chain-lighten + atom-lighten

        entity = array.entity_id[0]
        color_entity = self.color_entity[entity]
        nc = len(self.entity_chains[entity])
        ci = np.where(self.entity_chains[entity] == array.chain_name)[0][0] * 2
        color_chain = color.Lab.lighten_color(color_entity, (float(ci) / nc))
        colors = np.tile(color_chain, (len(array), 1))

        store_named_attribute(
            obj=obj,
            name="Color",
            data=colors,
            atype=AttributeTypes.FLOAT_COLOR,
        )

    def _create_object_instances(
        self, name: str = "CellPack", node_setup: bool = True
    ) -> bpy.types.Collection:
        collection = bl.coll.cellpack(name)

        for i, mol_id in enumerate(self.file.mol_ids):
            array = self.molecules[mol_id]

            obj = create_object(
                array=array,
                name=mol_id,
                collection=collection,
            )

            if len(self.color_entity) > 0:
                self._assign_colors(obj, array)

            if node_setup:
                nodes.create_starting_node_tree(
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
        mod = nodes.get_mod(self.object)

        group = nodes.new_tree(name=f"MN_ensemble_{name}", fallback=False)
        mod.node_group = group

        node_pack = nodes.add_custom(group, "Ensemble Instance", location=[-100, 0])
        node_pack.inputs["Instances"].default_value = self.data_collection
        node_pack.inputs["Fraction"].default_value = fraction
        node_pack.inputs["As Points"].default_value = as_points

        link = group.links.new
        link(nodes.get_input(group).outputs[0], node_pack.inputs[0])
        link(node_pack.outputs[0], nodes.get_output(group).inputs[0])
