from abc import ABCMeta
import bpy
from enum import Enum
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
    INTERACTION = "interaction"


class MolecularEntity(
    BlenderObject,
    metaclass=ABCMeta,
):
    def __init__(self) -> None:
        super().__init__(obj=None)
        self._entity_type: EntityType
        self._register_with_session()

    @property
    def bob(self) -> BlenderObject:
        return BlenderObject(self.object)

    @property
    def node_group(self) -> bpy.types.NodeGroup:
        return self.object.modifiers["MolecularNodes"].node_group

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
