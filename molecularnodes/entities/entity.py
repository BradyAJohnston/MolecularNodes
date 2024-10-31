from abc import ABCMeta
import bpy
from uuid import uuid1
from .. import blender as bl
from ..blender import databpy as db
from ..blender.databpy import AttributeTypes, AttributeType, Domains, DomainType
import warnings
import numpy as np


class ObjectMissingError(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class MolecularEntity(metaclass=ABCMeta):
    def __init__(self) -> None:
        self.uuid: str = str(uuid1())
        self._object: bpy.types.Object | None
        self.type: str = ""

    @property
    def name(self) -> str:
        obj = self.object
        if obj is None:
            return None

        return obj.name

    @name.setter
    def name(self, value: str) -> None:
        obj = self.object
        if obj is None:
            raise ObjectMissingError
        obj.name = value

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

    def named_attribute(self, name="position", evaluate=False) -> np.ndarray | None:
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
                "No object yet created. Use `create_object()` to create a corresponding object."
            )
            return None
        return db.named_attribute(self.object, name=name, evaluate=evaluate)

    def set_position(self, positions: np.ndarray) -> None:
        "A slightly optimised way to set the positions of the object's mesh"
        self.store_named_attribute(
            data=positions, name="position", atype=AttributeTypes.FLOAT_VECTOR
        )

    def set_boolean(self, boolean: np.ndarray, name="boolean") -> None:
        self.store_named_attribute(boolean, name=name, atype=AttributeTypes.BOOLEAN)

    def store_named_attribute(
        self,
        data: np.ndarray,
        name: str = "NewAttribute",
        atype: str | AttributeType | None = None,
        domain: str | DomainType = Domains.POINT,
        overwrite: bool = True,
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
                "No object yet created. Use `create_object()` to create a corresponding object."
            )
            return None
        db.store_named_attribute(
            obj=self.object,
            data=data,
            name=name,
            atype=atype,
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
            return list(bl.mesh.evaluate_object(cls.object).data.attributes.keys())

        return list(cls.object.data.attributes.keys())
