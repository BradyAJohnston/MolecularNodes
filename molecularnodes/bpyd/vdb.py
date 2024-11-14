import bpy
from .object import ObjectTracker
from .collection import create_collection


def import_vdb(
    file: str, collection: str | bpy.types.Collection | None = None
) -> bpy.types.Object:
    """
    Imports a VDB file as a Blender volume object.

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
        # Move the object to the given collection
        initial_collection = obj.users_collection[0]
        initial_collection.objects.unlink(obj)
        if isinstance(collection, str):
            collection = create_collection(collection)
        collection.objects.link(obj)

    return obj
