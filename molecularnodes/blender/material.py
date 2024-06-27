import bpy
import os
from .nodes import MN_DATA_FILE

materials = ["MN Default", "MN Flat Outline", "MN Squishy", "MN Transparent Outline"]


def append_material(name: str) -> bpy.types.Material:
    mat = bpy.data.materials.get(name)

    if not mat:
        bpy.ops.wm.append(
            directory=os.path.join(MN_DATA_FILE, "Material"),
            filename=name,
            link=False,
        )

    return bpy.data.materials[name]


def add_all_materials() -> None:
    for mat in materials:
        append_material(mat)


def default() -> bpy.types.Material:
    return append_material("MN Default")
