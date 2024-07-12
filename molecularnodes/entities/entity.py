from abc import ABCMeta
import bpy
from uuid import uuid1
from .. import blender as bl
import warnings
import numpy as np


class ObjectMissingError(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class MolecularEntity(metaclass=ABCMeta):
    def __init__(self) -> None:
        self.uuid: str = str(uuid1())
        self.object_ref: bpy.types.Object | None
        self.type: str = ""

    @property
    def name(self) -> str:
        bob = self.object
        if bob is None:
            return None

        return bob.name

    @name.setter
    def name(self, value: str) -> None:
        bob = self.object
        if bob is None:
            raise ObjectMissingError
        bob.name = value

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

    def get_attribute(self, name="position", evaluate=False) -> np.ndarray | None:
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
        if self.object is None:
            warnings.warn(
                "No object yet created. Use `create_model()` to create a corresponding object."
            )
            return None
        return bl.mesh.get_attribute(self.object, name=name, evaluate=evaluate)

    def set_position(self, positions: np.ndarray) -> None:
        "A slightly optimised way to set the positions of the object's mesh"
        bob = self.object
        attribute = bob.data.attributes["position"]
        n_points = len(attribute.data)
        if positions.shape != (n_points, 3):
            raise AttributeError(
                f"Expected an array of dimension {(n_points, 3)} to set the position"
                / f"but got {positions.shape=}"
            )

        # actually set the data for the positions
        attribute.data.foreach_set("vector", positions.reshape(-1))
        # trigger a depsgraph update. The second method is better but bugs out sometimes
        # so we try the first method initially
        try:
            bob.data.vertices[0].co = bob.data.vertices[0].co  # type: ignore
        except AttributeError:
            bob.data.update()  # type: ignore

    def set_boolean(self, boolean: np.ndarray, name="boolean") -> None:
        self.set_attribute(boolean, name=name, data_type="BOOLEAN")

    def set_attribute(
        self,
        data: np.ndarray,
        name="NewAttribute",
        data_type=None,
        domain="POINT",
        overwrite=True,
    ):
        """
        Set an attribute for the molecule.

        Parameters
        ----------
        data : np.ndarray
            The data to be set as the attribute. Must be of length equal to the length
            of the domain.
        name : str, optional
            The name of the new attribute. Default is 'NewAttribute'.
        type : str, optional
            If value is None (Default), the data type is inferred. The data type of the
            attribute. Possbible values are ('FLOAT_VECTOR', 'FLOAT_COLOR", 'QUATERNION',
            'FLOAT', 'INT', 'BOOLEAN').
        domain : str, optional
            The domain of the attribute. Default is 'POINT'. Possible values are
            currently ['POINT', 'EDGE', 'FACE', 'SPLINE']
        overwrite : bool, optional
            Whether to overwrite an existing attribute with the same name, or create a
            new attribute with always a unique name. Default is True.
        """
        if not self.object:
            warnings.warn(
                "No object yet created. Use `create_model()` to create a corresponding object."
            )
            return None
        bl.mesh.set_attribute(
            self.object,
            name=name,
            data=data,
            data_type=data_type,
            domain=domain,
            overwrite=overwrite,
        )

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
            return list(bl.mesh.evaluated(cls.object).data.attributes.keys())

        return list(cls.object.data.attributes.keys())
