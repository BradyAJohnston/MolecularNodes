import bpy
import databpy as db
import numpy as np
from . import coll


def create_data_object(
    array: np.ndarray,
    name: str = "DataObject",
    collection: bpy.types.Collection | None = None,
    world_scale: float = 0.1,
) -> bpy.types.Object:
    bob = db.create_bob(
        array["transform"][:, :3, 3] * world_scale,
        collection=collection if collection is not None else coll.data(),
        name=name,
    )

    attributes = [
        ("transform", db.mat4),
        ("assembly_id", db.int32),
        ("sym_id", db.int32),
        ("chain_id", db.int32),
        ("pdb_model_num", db.int32),
    ]

    for column, atype in attributes:
        try:
            data = array[column]
            if column == "transform":
                data[:, :3, 3] *= world_scale
        except ValueError:
            continue
        # us the unique sorted integer encoding version of the non-numeric
        # attribute, as GN doesn't support strings currently
        if np.issubdtype(data.dtype, str):
            data = np.unique(data, return_inverse=True)[1]

        bob[column] = atype(data)

    return bob.object
