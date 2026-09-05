from pathlib import Path
from typing import cast
import bpy
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
    ) -> bpy.types.Object:
        self._create_data_object(name=name)
        collection = self._create_object_instances(name=name, node_setup=node_setup)

        with self.tree.reset() as (atoms, join):
            (atoms >> geometry.EnsembleInstance(instances=collection) >> join)

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
            instance_tree.tree.is_modifier = True

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

        # stored by name via the `Ensemble` property, not as a live `Collection`
        # reference, which the session could not pickle
        self.instance_collection = collection

        return collection

    def _create_data_object(self, name: str = "DataObject") -> None:
        self.object = bl.mesh.create_data_object(
            self.transformations, name=name, collection=bl.coll.mn()
        )
        self.props.chain_ids = self.file.mol_ids.tolist()
        self.props.entity_type = self._entity_type.value

        # if we are dealing with petworld data, overwrite the chain_id for the data object
        if self.file._is_petworld:
            self.store_named_attribute(
                self.named_attribute("pdb_model_num"), "chain_id"
            )
