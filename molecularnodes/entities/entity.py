from abc import ABCMeta
import bpy
from uuid import uuid1
from .. import blender as bl
from ..blender.bpyd import (
    AttributeTypes,
    BlenderObject,
)
import warnings
import numpy as np


class MolecularEntity(
    BlenderObject,
    metaclass=ABCMeta,
):
    def __init__(self) -> None:
        self.uuid: str = str(uuid1())
        self._object: bpy.types.Object | None
        self.type: str = ""

    @property
    def object(self) -> bpy.types.Object | None:
        # If we don't have connection to an object, attempt to re-stablish to a new
        # object in the scene with the same UUID. This helps if duplicating / deleting
        # objects in the scene, but sometimes Blender just loses reference to the object
        # we are working with because we are manually setting the data on the mesh,
        # which can wreak havoc on the object database. To protect against this,
        # if we have a broken link we just attempt to find a new suitable object for it
        try:
            # if the connection is broken then trying to the name will raise a connection
            # error. If we are loading from a saved session then the object_ref will be
            # None and get an AttributeError
            self._object.name
            return self._object
        except (ReferenceError, AttributeError):
            for obj in bpy.data.objects:
                if obj.mn.uuid == self.uuid:
                    print(
                        Warning(
                            f"Lost connection to object: {self._object}, now connected to {obj}"
                        )
                    )
                    self._object = obj
                    return obj

            return None

    @object.setter
    def object(self, value):
        if isinstance(value, bpy.types.Object) or value is None:
            self._object = value
        else:
            raise TypeError(f"The `object` must be a Blender object, not {value=}")

    @property
    def bob(self) -> BlenderObject:
        return BlenderObject(self.object)

    def set_position(self, positions: np.ndarray) -> None:
        "A slightly optimised way to set the positions of the object's mesh"
        self.store_named_attribute(
            data=positions, name="position", atype=AttributeTypes.FLOAT_VECTOR
        )

    def set_boolean(self, boolean: np.ndarray, name="boolean") -> None:
        self.store_named_attribute(boolean, name=name, atype=AttributeTypes.BOOLEAN)

    @classmethod
    def list_attributes(cls, evaluate=False) -> list | None:
        """
        Returns a list of attribute names for the object.

        Parameters
        ----------
        evaluate : bool, optional
            Whether to first evaluate the modifiers on the object before listing the
            available attributes.

        Returns
        -------
        list[str] | None
            A list of attribute names if the molecule object exists, None otherwise.
        """
        if not cls.object:
            warnings.warn("No object created")
            return None
        if evaluate:
            return list(bl.mesh.evaluate_object(cls.object).data.attributes.keys())

        return list(cls.object.data.attributes.keys())
