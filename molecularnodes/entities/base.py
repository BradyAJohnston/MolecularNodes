from abc import ABCMeta
from enum import Enum
from typing import List
import bpy
from bpy.types import GeometryNodeTree
from databpy import (
    BlenderObject,
)
from nodebpy import geometry as g
from nodebpy.builder import GeometrySocket, TreeBuilder
from ..blender import utils as blender_utils
from ..nodes import nodes
from ..nodes.geometry import (
    GeometryNodeInterFace,
    style_interfaces_from_tree,
)
from .utilities import BoolObjectMNProperty


# create a EntityType enum
# These should match the entity_type EnumProperty in props.py
class EntityType(Enum):
    MD = "md"
    MD_OXDNA = "md-oxdna"
    MD_STREAMING = "md-streaming"
    MOLECULE = "molecule"
    ENSEMBLE = "ensemble"
    DENSITY = "density"
    ENSEMBLE_STAR = "ensemble-star"
    ENSEMBLE_CELLPACK = "ensemble-cellpack"


class MolecularTree(TreeBuilder[GeometryNodeTree]):
    def __init__(self, entitiy: "MolecularEntity", tree: TreeBuilder | str) -> None:
        self._entity = entitiy
        super().__init__(tree)

    def reset(self) -> tuple[GeometrySocket, GeometrySocket]:
        """Reset the tree to a default state. Returns the geometry input and the a JoinGeometry.i.geometry (input, join)"""
        self.clear()
        geo = self.inputs.geometry("Atoms")
        join = g.JoinGeometry()
        join >> self.outputs.geometry()
        return geo, join.i.geometry

    def clear(self) -> None:
        self.tree.nodes.clear()
        assert self.tree.interface
        self.tree.interface.clear()


class MolecularEntity(
    BlenderObject,
    metaclass=ABCMeta,
):
    update_with_scene = BoolObjectMNProperty("update_with_scene")
    _entity_type: EntityType | None = None

    def __init__(self) -> None:
        super().__init__(obj=None)
        self._register_with_session()
        self._world_scale = 0.01
        self._tree = None

    @property
    def node_group(self) -> bpy.types.GeometryNodeTree:
        mod = self.modifier
        mod.node_group = self.tree.tree
        return mod.node_group

    @property
    def modifier(self) -> bpy.types.NodesModifier:
        for mod in self.object.modifiers:
            if mod.type == "NODES" and mod.name == "Molecular Nodes":
                return mod

        mod = self.object.modifiers.new("GeometryNodes", "NODES")
        return mod

    @property
    def tree(self) -> MolecularTree:
        if self._tree is None:
            self._tree = MolecularTree(self, self.modifier_node_tree)
        return self._tree

    @property
    def modifier_node_tree(self) -> bpy.types.GeometryNodeTree:
        mod: bpy.types.NodesModifier = self.object.modifiers.get("Molecular Nodes")
        if mod is None:
            mod = self.object.modifiers.new("Molecular Nodes", "NODES")
        if mod.node_group is None:
            mod.node_group = g.tree().tree
            # raise ValueError(
            #     f"Unable to get MolecularNodes modifier for {self.object}, modifiers: {list(self.object.modifiers)}"
            # )
        return mod.node_group  # type: ignore

    @property
    def styles(self) -> List[GeometryNodeInterFace]:
        """
        Get the styles in the tree.
        """
        return style_interfaces_from_tree(self.modifier_node_tree)

    def _register_with_session(self) -> None:
        bpy.context.scene.MNSession.register_entity(self)  # type: ignore

    def set_frame(self, frame: int) -> None:
        """
        Update the underlying data to correspond to changes in the scene frame.

        This method should be implemented by subclasses to update the entity's
        data based on the given frame number. This can include updating positions
        and performing other necessary calculations.

        Args:
            frame (int): The frame number to update the entity's data to.

        Raises:
            NotImplementedError: If the method is not implemented by a subclass.
        """
        raise NotImplementedError("Subclasses must implement this method")

    def _setup_modifiers(self):
        """
        Create the modifiers for the molecule.
        """
        self.object.modifiers.new("Molecular Nodes", "NODES")
        # fallback=False => new tree all the time
        tree = nodes.new_tree(
            name=f"MN_{self.name}", input_name="Atoms", is_modifier=True, fallback=False
        )
        self.object.modifiers[0].node_group = tree

    def get_view(self) -> List[tuple]:
        """
        Get the 3D bounding box of the entity object

        """
        return blender_utils.get_bounding_box(self.object)

    def _get_annotation_entity_type(self) -> str:
        """
        Internal: Get entity type string for annotations.

        Subentities that derive from other entities can set this
        to the parent entity type to re-use annotations of the parent.
        By default, this returns the string value of the Entity type.

        Eg: OXDNA and StreamingTrajectory that derive from Trajectory entity
        can return EntityType.MD.value to re-use the annotations from
        the parent Trajectory entity.

        """
        return self._entity_type.value
