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
        self.type: str = ""
        self._object: bpy.types.Object | None

    @property
    def bob(self) -> BlenderObject:
        return BlenderObject(self.object)
