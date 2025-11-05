from abc import ABCMeta
from enum import Enum
from typing import List
import bpy
from databpy import (
    BlenderObject,
)
from ..blender import utils as blender_utils
from ..nodes import nodes
from ..nodes.geometry import (
    GeometryNodeInterFace,
    style_interfaces_from_tree,
)
from .utilities import BoolObjectMNProperty, _get_gn_modifier


# create a EntityType enum
# These should match the entity_type EnumProperty in props.py
class EntityType(Enum):
    MD = "md"
    MD_OXDNA = "md-oxdna"
    MOLECULE = "molecule"
    ENSEMBLE = "ensemble"
    DENSITY = "density"
    ENSEMBLE_STAR = "ensemble-star"
    ENSEMBLE_CELLPACK = "ensemble-cellpack"


bpy.types.NodesModifier


class MolecularEntity(
    BlenderObject,
    metaclass=ABCMeta,
):
    update_with_scene = BoolObjectMNProperty("update_with_scene")

    def __init__(self) -> None:
        super().__init__(obj=None)
        self._register_with_session()
        self._world_scale = 0.01

    @property
    def node_group(self) -> bpy.types.NodeTree | None:
        if "MolecularNodes" in self.object.modifiers:
            return _get_gn_modifier(self.object, "MolecularNodes").node_group
        return None

    @property
    def tree(self) -> bpy.types.GeometryNodeTree:
        mod: bpy.types.NodesModifier = self.object.modifiers["MolecularNodes"]  # type: ignore
        if mod is None:
            raise ValueError(
                f"Unable to get MolecularNodes modifier for {self.object}, modifiers: {list(self.object.modifiers)}"
            )
        return mod.node_group  # type: ignore

    @property
    def styles(self) -> List[GeometryNodeInterFace]:
        """
        Get the styles in the tree.
        """
        return style_interfaces_from_tree(self.tree)

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
        self.object.modifiers.new("MolecularNodes", "NODES")
        # fallback=False => new tree all the time
        tree = nodes.new_tree(  # type: ignore
            name=f"MN_{self.name}", input_name="Atoms", is_modifier=True, fallback=False
        )
        self.object.modifiers[0].node_group = tree  # type: ignore

    def get_view(self) -> List[tuple]:
        """
        Get the 3D bounding box of the entity object

        """
        return blender_utils.get_bounding_box(self.object)
