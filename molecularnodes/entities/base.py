from abc import ABCMeta
from enum import Enum
import bpy
from databpy import (
    BlenderObject,
)


# create a EntityType enum with strings for values, "md", "md-oxdna", "molecule", "star"
class EntityType(Enum):
    MD = "md"
    MD_OXDNA = "md-oxdna"
    MOLECULE = "molecule"
    ENSEMBLE = "ensemble"
    DENSITY = "density"


class MolecularEntity(
    BlenderObject,
    metaclass=ABCMeta,
):
    def __init__(self) -> None:
        super().__init__(obj=None)
        self._entity_type: EntityType
        self._register_with_session()
        self._world_scale = 0.01

    @property
    def bob(self) -> BlenderObject:
        return BlenderObject(self.object)

    @property
    def node_group(self) -> bpy.types.GeometryNodeTree:
        return self.object.modifiers["MolecularNodes"].node_group

    @property
    def update_with_scene(self) -> bool:
        return self.object.mn.update_with_scene

    @update_with_scene.setter
    def update_with_scene(self, value: bool) -> None:
        self.object.mn.update_with_scene = value

    def _register_with_session(self) -> None:
        bpy.context.scene.MNSession.register_entity(self)

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
