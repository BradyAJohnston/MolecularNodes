from pathlib import Path
from typing import cast
import bpy
from databpy import BlenderObject
from nodebpy.builder import TreeBuilder
from ... import blender as bl
from ...nodes import geometry, material
from ..utilities import create_object
from .base import Ensemble, EntityType
from .reader import CellPackReader


class CellPack(Ensemble):
    def __init__(self, file_path: str | Path) -> None:
        super().__init__(file_path)
        self._entity_type = EntityType.ENSEMBLE_CELLPACK
        self.file = CellPackReader(file_path)
        self.transformations = self.file.get_assemblies()

    @property
    def molecules(self) -> dict:
        return self.file.molecules

    def create_object(
        self,
        name: str = "CellPack",
        node_setup: bool = True,
        world_scale: float = 0.01,
        fraction: float = 1.0,
        simplify: bool = False,
    ) -> bpy.types.Object:
        self.object = self._create_data_object(name=name)
        self.object.mn.entity_type = self._entity_type.value
        self._create_object_instances(name=name, node_setup=node_setup)

        with self.tree.reset() as (atoms, join):
            (
                atoms
                >> geometry.EnsembleInstance(
                    instances=self.data_collection, fraction=fraction
                )
                >> join
            )

        return self.object

    def _create_object_instances(
        self, name: str = "CellPack", node_setup: bool = True
    ) -> bpy.types.Collection:
        collection = bl.coll.cellpack(name)

        # a single node group is shared by every instanced molecule rather than
        # building an identical tree per-object
        if node_setup:
            with TreeBuilder.geometry(f"MN_pack_instance_{name}") as instance_tree:
                atoms = instance_tree.inputs.geometry("Atoms")
                style = geometry.StyleSpheres()
                material.assign_material(style.node, "MN Ambient Occlusion")
                atoms >> style >> instance_tree.outputs.geometry("Geometry")

        for mol_id in self.file.mol_ids:
            obj = create_object(
                array=self.molecules[mol_id],
                name=mol_id,
                collection=collection,
            )
            obj.mn.entity_type = self._entity_type.value

            if node_setup:
                mod = cast(
                    bpy.types.NodesModifier,
                    obj.modifiers.new("Molecular Nodes", "NODES"),
                )
                mod.node_group = instance_tree.tree

        self.data_collection = collection
        self.instance_collection = collection

        return collection

    def _create_data_object(self, name: str = "DataObject") -> bpy.types.Object:
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
