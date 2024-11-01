from abc import ABCMeta
import bpy
from uuid import uuid1
from .. import blender as bl
from ..blender.bpyd import (
    AttributeTypes,
    BlenderObject,
)
import numpy as np


class MolecularEntity(
    BlenderObject,
    metaclass=ABCMeta,
):
    def __init__(self) -> None:
        self.uuid: str = str(uuid1())
        self._object: bpy.types.Object | None
        self.type: str = ""

    @object.setter
    def object(self, value):
        if isinstance(value, bpy.types.Object) or value is None:
            self._object = value
        else:
            raise TypeError(f"The `object` must be a Blender object, not {value=}")

    @property
    def bob(self) -> BlenderObject:
        return BlenderObject(self.object)
