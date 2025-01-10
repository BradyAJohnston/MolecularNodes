from bpy.types import Material
from databpy.material import append_from_blend
from ..utils import MN_DATA_FILE

MATERIAL_NAMES = [
    "MN Default",
    "MN Flat Outline",
    "MN Squishy",
    "MN Transparent Outline",
    "MN Ambient Occlusion",
]


def append(name: str) -> Material:
    "Append a material from the MN_DATA_FILE."
    return append_from_blend(name, MN_DATA_FILE)


def add_all_materials() -> None:
    "Append all pre-defined materials from the MN_DATA_FILE."
    for name in MATERIAL_NAMES:
        append(name)


def default() -> Material:
    "Return the default material."
    return append("MN Default")
