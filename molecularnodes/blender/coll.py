import bpy
from typing import Optional


def mn() -> bpy.types.Collection:
    """Return the MolecularNodes Collection

    The collection called 'MolecularNodes' inside the Blender scene is returned. If the
    collection does not exist first, it is created.
    """
    coll = bpy.data.collections.get("MolecularNodes")
    if not coll:
        coll = bpy.data.collections.new("MolecularNodes")
        bpy.context.scene.collection.children.link(coll)
    return coll


def armature(name: str = "MN_armature") -> bpy.types.Collection:
    coll = bpy.data.collections.get(name)
    if not coll:
        coll = bpy.data.collections.new(name)
        mn().children.link(coll)
    return coll


def data(suffix: str = "") -> bpy.types.Collection:
    """A collection for storing MN related data objects."""
    name = f"MN_data{suffix}"

    collection = bpy.data.collections.get(name)
    if not collection:
        collection = bpy.data.collections.new(name)
        mn().children.link(collection)

        # disable the view of the data collection
        bpy.context.view_layer.layer_collection.children["MolecularNodes"].children[name].exclude = True
    return collection


def frames(name: str = "", parent: Optional[bpy.types.Object] = None, suffix: str = "_frames") -> bpy.types.Collection:
    """Create a Collection for Frames of a Trajectory

    Args:
        name (str, optional): Name of the collection for the frames. Defaults to "".
        parent (_type_, optional): A blender collection which will become the parent
        collection. Defaults to the MolecularNodes collection if None.
    """
    coll_frames = bpy.data.collections.new(name + suffix)
    if not parent:
        mn().children.link(coll_frames)
    else:
        parent.children.link(coll_frames)

    return coll_frames


def cellpack(name: str = "", parent: Optional[bpy.types.Object] = None, fallback: bool = False) -> bpy.types.Collection:
    full_name = f"cellpack_{name}"
    coll = bpy.data.collections.get(full_name)
    if coll and fallback:
        return coll

    coll = bpy.data.collections.new(full_name)

    if parent:
        parent.children.link(coll)
    else:
        data().children.link(coll)

    return coll
