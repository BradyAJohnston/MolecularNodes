from abc import ABCMeta
import bpy
from uuid import uuid1
from . import blender as bl
import warnings
import numpy as np


class ObjectMissingError(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class MNDataObject(metaclass=ABCMeta):
    def __init__(self) -> None:
        self.name: str | None
        self.object = None
        self.object_ref = None
        self.uuid: str = str(uuid1())

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
            self.object_ref.name
            return self.object_ref
        except (ReferenceError, AttributeError):
            for bob in bpy.data.objects:
                if bob.mn.uuid == self.uuid:
                    print(
                        Warning(
                            f"Lost connection to object: {self.object_ref}, now connected to {bob}"
                        )
                    )
                    self.object_ref = bob
                    return bob

            return None

    @object.setter
    def object(self, value):
        self.object_ref = value

    @classmethod
    def get_attribute(cls, name="position", evaluate=False) -> np.ndarray | None:
        """
        Get the value of an object for the data molecule.

        Parameters
        ----------
        name : str, optional
            The name of the attribute. Default is 'position'.
        evaluate : bool, optional
            Whether to first evaluate all node trees before getting the requsted attribute.
            False (default) will sample the underlying atomic geometry, while True will
            sample the geometry that is created through the Geometry Nodes tree.

        Returns
        -------
        np.ndarray
            The value of the attribute.
        """
        if not cls.object:
            warnings.warn(
                "No object yet created. Use `create_model()` to create a corresponding object."
            )
            return None
        return bl.obj.get_attribute(cls.object, name=name, evaluate=evaluate)

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
            return list(bl.obj.evaluated(cls.object).data.attributes.keys())

        return list(cls.object.data.attributes.keys())
