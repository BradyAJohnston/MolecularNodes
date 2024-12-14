from abc import ABCMeta
import bpy
from enum import Enum
from uuid import uuid1
from ..bpyd import (
    BlenderObject,
)


# create a EntityType enum with strings for values, "md", "md-oxdna", "molecule", "star"
class EntityType(Enum):
    MD = "md"
    MD_OXDNA = "md-oxdna"
    MOLECULE = "molecule"
    STAR = "star"


class MolecularEntity(
    BlenderObject,
    metaclass=ABCMeta,
):
    def __init__(self) -> None:
        self.uuid: str = str(uuid1())
        self.type: str = ""
        self._object: bpy.types.Object | None

    @property
    def bob(self) -> BlenderObject:
        return BlenderObject(self.object)

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
