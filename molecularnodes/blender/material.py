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


# class to interact with a bpy.types.Material node tree and change some of the default
# values of the nodes inside of it
class MaterialTreeInterface:
    def __init__(self, tree: Material):
        self.tree = tree
        self.nodes = tree.node_tree.nodes
        self.links = tree.node_tree.links


class MaterialAmbientOcclusion(MaterialTreeInterface):
    def __init__(self, tree: Material):
        super().__init__(tree)

    @property
    def ao_power(self) -> float:
        return self.nodes["Math"].inputs[1].default_value

    @ao_power.setter
    def ao_power(self, value: float):
        self.nodes["Math"].inputs[1].default_value = value

    @property
    def ao_distance(self) -> float:
        return self.nodes["Ambient Occlusion"].inputs["Distance"].default_value

    @ao_distance.setter
    def ao_distance(self, value: float) -> None:
        self.nodes["Ambient Occlusion"].inputs["Distance"].default_value = value


class MaterialSquishy(MaterialTreeInterface):
    def __init__(self, tree: Material):
        super().__init__(tree)

    @property
    def subsurf_scale(self) -> float:
        return self.nodes["Principled BSDF"].inputs["Subsurface Scale"].default_value

    @subsurf_scale.setter
    def subsurf_scale(self, value: float):
        self.nodes["Principled BSDF"].inputs["Subsurface Scale"].default_value = value
