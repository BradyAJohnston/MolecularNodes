from bpy.types import Material
from databpy.material import append_from_blend
from ..assets import MN_DATA_FILE
from typing import Any, TypeVar, Generic, Union

MATERIAL_NAMES = [
    "MN Default",
    "MN Flat Outline",
    "MN Squishy",
    "MN Transparent Outline",
    "MN Ambient Occlusion",
]


def append(name: str) -> Material:
    "Append a material from the MN_DATA_FILE."
    return append_from_blend(name, str(MN_DATA_FILE))


def add_all_materials() -> None:
    "Append all pre-defined materials from the MN_DATA_FILE."
    for name in MATERIAL_NAMES:
        append(name)


def default() -> Material:
    "Return the default material."
    return append("MN Default")


# class to interact with a bpy.types.Material node tree and change some of the default
# values of the nodes inside of it
class MaterialTreeInterface:
    def __init__(self):
        self._material: Material

    @property
    def material(self):
        return self._material

    @property
    def node_tree(self):
        return self._material.node_tree

    @property
    def nodes(self):
        return self.node_tree.nodes

    @property
    def links(self):
        return self.node_tree.links


T = TypeVar("T", float, int, tuple[float, ...], bool)


def socket(node_name: str, socket: str | int, type_: type[T]):
    def getter(self) -> T:
        return self.nodes[node_name].inputs[socket].default_value

    def setter(self, value: T) -> None:
        self.nodes[node_name].inputs[socket].default_value = value

    return property(getter, setter)


def option(node_name: str, input: str | int, type_: type[T]):
    def getter(self) -> T:
        return self.nodes[node_name][input].default_value

    def setter(self, value: T) -> None:
        self.nodes[node_name][input].default_value = value

    return property(getter, setter)


class AmbientOcclusion(MaterialTreeInterface):
    def __init__(self):
        self._material = append("MN Ambient Occlusion")

    power = socket("Math", 1, float)
    distance = socket("Ambient Occlusion", "Distance", float)
    samples = socket("Ambient Occlusion", "Samples", int)


class MaterialSquishy(MaterialTreeInterface):
    def __init__(self):
        self._material = append("MN Squishy")

    subsurface_scale = socket("Principled BSDF", "Subsurface Scale", float)


class Default(MaterialTreeInterface):
    def __init__(self):
        self._material = append("MN Default")

    metallic = socket("Principled BSDF", "Metallic", float)
    roughness = socket("Principled BSDF", "Roughness", float)
    ior = socket("Principled BSDF", "IOR", float)
    alpha = socket("Principled BSDF", "Alpha", float)

    diffuse_roughness = socket("Principled BSDF", "Diffuse Roughness", float)
    subsurface_weight = socket("Principled BSDF", "Subsurface Weight", float)
    subsurface_radius = socket(
        "Principled BSDF", "Subsurface Radius", tuple[float, float, float]
    )
    subsurface_scale = socket("Principled BSDF", "Subsurface Scale", float)
