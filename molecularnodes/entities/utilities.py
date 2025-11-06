"""
Python module for converting molecular structures from biotite into Blender objects.

This module provides functions to:
1. Create 3D models from biotite AtomArray objects (`create_object`)
2. Attach molecular information as attributes to Blender objects (`atom_array_to_named_attributes`)

The created Blender objects can be used for visualization and manipulation of molecular
structures within Blender's 3D environment. Properties like atom positions, bonds, and
other molecular attributes are preserved as named attributes that can be used with
Geometry Nodes and other Blender features.
"""

from typing import Any, Callable
import bpy
import databpy as db
import numpy as np
from biotite.structure import AtomArray, AtomArrayStack

__all__ = ["create_object", "atom_array_to_named_attributes"]


def _get_gn_modifier(obj: bpy.types.Object, name: str) -> bpy.types.NodesModifier:
    mod = obj.modifiers[name]
    if not isinstance(mod, bpy.types.NodesModifier):
        raise TypeError(
            f"Modifier {name} is of type {type(mod)} not {bpy.types.NodesModifier}"
        )
    return mod


def _unique_aname(obj: bpy.types.Object, prefix: str = "sel") -> str:
    attributes = db.list_attributes(obj)
    counter = 0
    aname = "{}_{}".format(prefix, counter)
    while aname in attributes:
        counter += 1
        aname = "{}_{}".format(prefix, counter)

    return aname


def _validate_non_negative(value: int) -> None:
    """Validate non-negative integers.

    Parameters
    ----------
    value : int
        Value to validate

    Raises
    ------
    ValueError
        If value is negative
    """
    if value < 0:
        raise ValueError(f"Value must be non-negative, got {value}")


class ObjectMNProperty:
    """Descriptor for Blender object properties with optional validation.

    Provides clean interface for Blender custom properties with validation support.

    Examples
    --------
    >>> class MyClass:
    ...     frame = ObjectMNProperty("frame", validate_fn=lambda x: x >= 0)
    """

    def __init__(
        self,
        attr_name: str,
        validate_fn: Callable[[Any], None] | None = None,
    ):
        """Initialize the descriptor.

        Parameters
        ----------
        attr_name : str
            Name of the property on the Blender object
        validate_fn : callable, optional
            Validation function that raises on invalid values
        """
        self.attr_name = attr_name
        self.validate_fn = validate_fn

    def __set_name__(self, owner, name):
        """Store the Python attribute name."""
        self.name = name

    def __get__(self, obj, objtype=None):
        """Get property value from Blender object."""
        return getattr(obj.object.mn, self.attr_name)

    def __set__(self, obj, value):
        """Set property value on Blender object with validation."""
        if self.validate_fn is not None:
            self.validate_fn(value)
        setattr(obj.object.mn, self.attr_name, value)


class IntObjectMNProperty(ObjectMNProperty):
    def __get__(self, obj, objtype=None) -> int:
        return super().__get__(obj, objtype)

    def __set__(self, obj, value: int):
        return super().__set__(obj, value)


class BoolObjectMNProperty(ObjectMNProperty):
    def __get__(self, obj, objtype=None) -> bool:
        return super().__get__(obj, objtype)

    def __set__(self, obj, value: bool):
        return super().__set__(obj, value)


class StringObjectMNProperty(ObjectMNProperty):
    def __get__(self, obj, objtype=None) -> str:
        return super().__get__(obj, objtype)

    def __set__(self, obj, value: str):
        return super().__set__(obj, value)


def create_object(
    array: AtomArray | AtomArrayStack,
    name="NewObject",
    centre: str | None = None,
    world_scale: float = 0.01,
    collection: bpy.types.Collection | None = None,
) -> bpy.types.Object:
    """
    Create a 3D model of the molecule, with one vertex for each atom.

    This function converts a biotite AtomArray into a Blender object with vertices
    representing atoms and edges representing bonds. All atomic properties from the
    AtomArray are stored as named attributes on the Blender mesh, making them
    accessible for use with Geometry Nodes and other Blender features.

    Parameters
    ----------
    array : AtomArray
        The molecular structure to convert into a Blender object.
    name : str, optional
        Name of the created Blender object. Default is "NewObject".
    centre : str or None, optional
        If provided, centers the object according to the specified criteria.
        Default is None (no centering).
    world_scale : float, optional
        Scale factor to apply to atomic coordinates when creating the model.
        Default is 0.01, which converts Ångströms to a reasonable size in Blender.
    collection : bpy.types.Collection or None, optional
        The Blender collection to add the new object to. If None, the object
        is added to the active collection. Default is None.

    Returns
    -------
    bpy.types.Object
        The created 3D model, as an object in the 3D scene.
        The object contains vertices for atoms and edges for bonds,
        with all molecular properties stored as named attributes.
    """

    if isinstance(array, AtomArrayStack):
        array = array[0]

    bob = db.create_bob(
        vertices=array.coord * world_scale,
        edges=array.bonds.as_array()[:, :2] if array.bonds is not None else None,
        name=name,
        collection=collection,
    )
    # Add information about the bond types to the model on the edge domain
    # Bond types: 'ANY' = 0, 'SINGLE' = 1, 'DOUBLE' = 2, 'TRIPLE' = 3, 'QUADRUPLE' = 4
    # 'AROMATIC_SINGLE' = 5, 'AROMATIC_DOUBLE' = 6, 'AROMATIC_TRIPLE' = 7
    # https://www.biotite-python.org/apidoc/biotite.structure.BondType.html#biotite.structure.BondType
    if array.bonds:
        bob.store_named_attribute(
            array.bonds.as_array()[:, 2],
            "bond_type",
            domain=db.AttributeDomains.EDGE,
            atype=db.AttributeTypes.INT,
        )

    atom_array_to_named_attributes(array, bob.object, world_scale=world_scale)

    if centre is not None:
        bob.position -= bob.centroid(centre)

    return bob.object


def atom_array_to_named_attributes(
    array: AtomArray, obj: bpy.types.Object, world_scale: float = 0.01
) -> None:
    """
    Store all annotations from an AtomArray as named attributes on Blender vertex data.

    This function iterates through the annotations (properties) of a biotite AtomArray
    and stores them as named attributes on a Blender mesh object. These attributes can
    then be accessed and manipulated using Blender's Geometry Nodes system or through Python.

    Only numeric and boolean attributes are stored, as Blender's Geometry Nodes doesn't
    currently support string attributes. String attributes from the AtomArray should have
    been converted to numeric representations during the import process.

    Parameters
    ----------
    array : AtomArray
        The biotite AtomArray containing molecular data and annotations to be stored.
    obj : bpy.types.Object
        The Blender object (typically a mesh) on which to store the attributes.
    world_scale : float, optional
        Scale factor to apply to size-related values like van der Waals radii.
        Default is 0.01 to convert from Ångströms to Blender units.

    Notes
    -----
    - Coordinate data ('coord') is skipped as it's already stored as vertex positions.
    - 'hetero' flag is also skipped from the annotations.
    - Any string attributes that were converted to integer codes will have the '_int'
      suffix removed when stored on the mesh.
    """

    # don't need to add coordinates as those have been stored as `position` on the mesh
    annotations_to_skip = ["coord", "hetero"]

    for attr in array.get_annotation_categories():
        if attr in annotations_to_skip:
            continue

        data = array.get_annotation(attr)  # type: ignore

        if attr == "vdw_radii":
            data *= world_scale

        # geometry nodes doesn't support strings at the moment so we can only store the
        # numeric and boolean attributes on the mesh. All string attributes should have
        # already been converted to a corresponding numeric attribute during the
        # reader process

        if not (
            np.issubdtype(data.dtype, np.number) or np.issubdtype(data.dtype, np.bool_)
        ):
            continue

        # the integer versions of strings have been added as annotations that just
        # append `_int` onto the name so we don't overrite the original data but when
        # storing on the meshh we can just remove the `_int`
        db.store_named_attribute(
            obj=obj,
            data=data,
            name=attr.replace("_int", ""),
        )
