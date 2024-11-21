from abc import ABCMeta
import bpy
from uuid import uuid1
from bpy.types import Object
from ..bpyd import (
    BlenderObject,
)


class MolecularEntity(
    BlenderObject,
    metaclass=ABCMeta,
):
    def __init__(self) -> None:
        self.uuid: str = str(uuid1())
        bpy.context.scene.MNSession.entities[self.uuid] = self

    @property
    def bob(self) -> BlenderObject:
        return BlenderObject(self.object)

    @property
    def object(self) -> Object:
        try:
            return bpy.data.objects[self._object_name]
        except KeyError:
            # if we can't find a refernce to the object via a name, then we look via the
            # unique uuids that were assigned to the object and the entity to match up
            for obj in bpy.data.objects:
                if obj.mn.uuid == self.uuid:
                    self._object_name = obj.name
                    return obj

    @object.setter
    def object(self, value) -> None:
        if not isinstance(value, Object):
            raise ValueError(
                f"Can only set object to be of type bpy.types.Object, not {type(value)=}"
            )
        self._object_name = value.name
