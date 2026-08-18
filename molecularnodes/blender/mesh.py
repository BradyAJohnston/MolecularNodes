import bpy
import numpy as np
from databpy.attribute import AttributeTypes
from databpy.object import create_bob
from ..nodes import nodes
from . import coll


def evaluate_using_mesh(obj: bpy.types.Object) -> bpy.types.Object:
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
    bob = create_bob()
    mod = nodes.get_mod(bob.object)
    mod.node_group = nodes.create_debug_group()
    mod.node_group.nodes["Object Info"].inputs["Object"].default_value = obj

    # need to use 'evaluate' otherwise the modifiers won't be taken into account
    return bob.evaluate()


def create_data_object(
    array: np.ndarray,
    name: str = "DataObject",
    collection: bpy.types.Collection | None = None,
    world_scale: float = 0.1,
) -> bpy.types.Object:
    bob = create_bob(
        array["transform"][:, :3, 3] * world_scale,
        collection=collection if collection is not None else coll.data(),
        name=name,
    )

    attributes = [
        ("transform", AttributeTypes.FLOAT4X4),
        ("assembly_id", AttributeTypes.INT),
        ("sym_id", AttributeTypes.INT),
        ("chain_id", AttributeTypes.INT),
        ("pdb_model_num", AttributeTypes.INT),
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

        # Blender stores float4x4 attributes column-major, so the row-major
        # numpy matrices must be transposed or GN sees the inverse rotation
        if type == AttributeTypes.FLOAT4X4:
            data = data.transpose(0, 2, 1)

        bob.store_named_attribute(data=data, name=column, atype=type)

    return bob.object
