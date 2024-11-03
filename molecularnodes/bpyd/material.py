from bpy.types import Material
import bpy
import os


# TODO: use DuplicatePrevention when adding material node trees
def append_from_blend(name: str, filepath: str) -> Material:
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Given file not found: {filepath}")
    try:
        return bpy.data.materials[name]
    except KeyError:
        bpy.ops.wm.append(
            directory=os.path.join(filepath, "Material"),
            filename=name,
            link=False,
        )
        return bpy.data.materials[name]
