from typing import Optional

import bpy
import numpy as np

from . import coll, nodes
from .databpy.attribute import (
    Attribute,
    AttributeMismatchError,
    AttributeTypes,
    guess_atype_from_array,
    store_named_attribute,
    named_attribute,
)
from .databpy.object import ObjectTracker, create_object
from .databpy.utils import evaluate_object


def centre(position: np.ndarray):
    "Calculate the centroid of the vectors"
    return np.mean(position, axis=0)


def centre_weighted(position: np.ndarray, weight: np.ndarray):
    "Calculate the weighted centroid of the vectors"
    return np.sum(position * weight.reshape((-1, 1)), axis=0) / np.sum(weight)


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
        ("rotation", AttributeTypes.QUATERNION),
        ("assembly_id", AttributeTypes.INT),
        ("chain_id", AttributeTypes.INT),
        ("transform_id", AttributeTypes.INT),
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

        store_named_attribute(obj=obj, data=data, name=column, atype=type)

    return obj
