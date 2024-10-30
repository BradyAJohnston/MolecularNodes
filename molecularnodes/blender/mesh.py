from typing import Optional

import bpy
import numpy as np

from . import coll, nodes
from .attribute import (
    Attribute,
    AttributeMismatchError,
    AttributeTypes,
    guess_atype_from_array,
)
from .object import ObjectTracker, create_object, evaluate_object


def centre(array: np.ndarray):
    return np.mean(array, axis=0)

def centre_weighted(array: np.ndarray, weight: np.ndarray):
    return np.sum(array * weight.reshape((len(array), 1)), axis=0) / np.sum(weight)


def store_named_attribute(
    obj: bpy.types.Object,
    data: np.ndarray,
    name: str,
    data_type: Optional[str] = None,
    domain: str = "POINT",
    overwrite: bool = True,
) -> bpy.types.Attribute:
    """
    Adds and sets the values of an attribute on the object.

    Parameters
    ----------
    obj : bpy.types.Object
        The Blender object.
    name : str
        The name of the attribute.
    data : np.ndarray
        The attribute data as a numpy array.
    type : str, optional
        The data type of the attribute. Defaults to None. Possible values are:
        'FLOAT_VECTOR', 'FLOAT_COLOR', 'FLOAT4X4', 'QUATERNION', 'FLOAT', 'INT', 'BOOLEAN'
    domain : str, optional
        The domain of the attribute. Defaults to 'POINT'. Currently, only 'POINT', 'EDGE',
        and 'FACE' have been tested.
    overwrite : bool, optional
        Whether to overwrite an existing attribute with the same name. Defaults to True.

    Returns
    -------
    bpy.types.Attribute
        The added attribute.
    """

    if isinstance(data_type, str):
        try:
            atype = AttributeTypes[data_type].value
        except KeyError:
            raise ValueError(
                f"Given data type {data_type=} does not match any of the possible attribute types: {list(AttributeTypes)=}"
            )

    if data_type is None:
        atype = guess_atype_from_array(data)

    attribute = obj.data.attributes.get(name)  # type: ignore
    if not attribute or not overwrite:
        attribute = obj.data.attributes.new(name, atype.type_name, domain)

    if len(data) != len(attribute.data):
        raise AttributeMismatchError(
            f"Data length {len(data)}, dimensions {data.shape} does not equal the size of the target domain {domain}, len={len(attribute.data)=}"
        )

    # the 'foreach_set' requires a 1D array, regardless of the shape of the attribute
    # it also requires the order to be 'c' or blender might crash!!
    attribute.data.foreach_set(atype.value_name, data.reshape(-1))

    # The updating of data doesn't work 100% of the time (see:
    # https://projects.blender.org/blender/blender/issues/118507) so this resetting of a
    # single vertex is the current fix. Not great as I can see it breaking when we are
    # missing a vertex - but for now we shouldn't be dealing with any situations where this
    # is the case For now we will set a single vert to it's own position, which triggers a
    # proper refresh of the object data.
    try:
        obj.data.vertices[0].co = obj.data.vertices[0].co  # type: ignore
    except AttributeError:
        obj.data.update()  # type: ignore

    return attribute


def named_attribute(
    obj: bpy.types.Object, name="position", evaluate=False
) -> np.ndarray:
    """
    Get the named attribute data from the object, optionally evaluating modifiers first.

    Parameters:
        object (bpy.types.Object): The Blender object.
        name (str, optional): The name of the attribute. Defaults to 'position'.

    Returns:
        np.ndarray: The attribute data as a numpy array.
    """
    if evaluate:
        obj = evaluate_object(obj)
    verbose = False
    try:
        attr = Attribute(obj.data.attributes[name])
    except KeyError:
        message = f"The selected attribute '{name}' does not exist on the mesh."
        if verbose:
            message += f"Possible attributes are: {obj.data.attributes.keys()}"

        raise AttributeError(message)

    return attr.as_array()


def import_vdb(file: str, collection: bpy.types.Collection = None) -> bpy.types.Object:
    """
    Imports a VDB file as a Blender volume object, in the MolecularNodes collection.

    Parameters
    ----------
    file : str
        Path to the VDB file.

    Returns
    -------
    bpy.types.Object
        A Blender object containing the imported volume data.
    """

    # import the volume object
    with ObjectTracker() as o:
        bpy.ops.object.volume_import(filepath=file, files=[])
        obj = o.latest()

    if collection:
        # Move the object to the MolecularNodes collection
        initial_collection = obj.users_collection[0]
        initial_collection.objects.unlink(obj)
        collection = coll.mn()
        collection.objects.link(obj)

    return obj




def evaluate_using_mesh(obj):
    """
    Evaluate the object using a debug object. Some objects can't currently have their
    Geometry Node trees evaluated (such as volumes), so we source the geometry they create
    into a mesh object, which can be evaluated and tested.

    Parameters
    ----------
    object : bpy.types.Object
        The object to be evaluated.

    Returns
    -------
    bpy.types.Object

    Notes
    -----
    Intended for debugging only.
    """
    # create an empty mesh object. It's modifiers can be evaluated but some other
    # object types can't be currently through the API
    debug_obj = create_object()
    mod = nodes.get_mod(debug_obj)
    mod.node_group = nodes.create_debug_group()
    mod.node_group.nodes["Object Info"].inputs["Object"].default_value = obj

    # need to use 'evaluate' otherwise the modifiers won't be taken into account
    return evaluate_object(debug_obj)


def create_data_object(array, collection=None, name="DataObject", world_scale=0.01):
    # still requires a unique call TODO: figure out why
    # I think this has to do with the bcif instancing extraction
    array = np.unique(array)
    locations = array["translation"] * world_scale

    if not collection:
        collection = coll.data()

    obj = create_object(locations, collection=collection, name=name)

    attributes = [
        ("rotation", "QUATERNION"),
        ("assembly_id", "INT"),
        ("chain_id", "INT"),
        ("transform_id", "INT"),
    ]

    for column, type in attributes:
        try:
            data = array[column]
        except ValueError:
            continue
        # us the unique sorted integer encoding version of the non-numeric
        # attribute, as GN doesn't support strings currently
        if np.issubdtype(data.dtype, str):
            data = np.unique(data, return_inverse=True)[1]

        store_named_attribute(obj=obj, data=data, name=column, data_type=type)

    return obj
