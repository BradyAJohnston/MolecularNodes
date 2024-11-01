import bpy
from bpy.types import Collection


def create_collection(
    name: str = "NewCollection", parent: Collection | str | None = None
) -> Collection:
    if isinstance(parent, str):
        try:
            parent = bpy.data.collections[name]
        except KeyError:
            parent = bpy.data.collections.new(name)
            bpy.context.scene.collection.children.linke(parent)
    try:
        coll = bpy.data.collections[name]
    except KeyError:
        coll = bpy.data.collections.new(name)
        if parent is None:
            bpy.context.scene.collection.children.link(coll)
        else:
            parent.children.link(coll)

    return coll
