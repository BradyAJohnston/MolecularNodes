from bpy.types import Material, ShaderNodeTree
from .nodes import TreeInterface
from .utils import append_material, socket, option
from .utils import MATERIAL_NAMES


def default() -> Material:
    "Return the default material."
    return append_material("MN Default")


# class to interact with a bpy.types.Material node tree and change some of the default
# values of the nodes inside of it
class MaterialTreeInterface(TreeInterface):
    def __init__(self):
        self._material: Material

    @property
    def material(self) -> Material:
        return self._material

    @property
    def node_tree(self) -> ShaderNodeTree:
        if self._material.node_tree is None:
            raise ValueError("Material has no node tree")
        return self._material.node_tree

    @property
    def nodes(self):
        return self.node_tree.nodes

    @property
    def links(self):
        return self.node_tree.links


def get_material_interface(material: str | Material) -> MaterialTreeInterface:
    if isinstance(material, Material):
        material = material.name
    match material:
        case "MN Default":
            return Default()
        case "MN Flat Outline":
            return FlatOutline()
        case "MN Squishy":
            return Squishy()
        case "MN Transparent Outline":
            return TransparentOutline()
        case "MN Ambient Occlusion":
            return AmbientOcclusion()
        case _:
            raise ValueError(f"Unknown material name: {material}")


class AmbientOcclusion(MaterialTreeInterface):
    def __init__(self):
        self._material = append_material("MN Ambient Occlusion")

    power = socket("Math", 1, float)
    distance = socket("Ambient Occlusion", "Distance", float)
    samples = socket("Ambient Occlusion", "Samples", int)


class FlatOutline(MaterialTreeInterface):
    def __init__(self):
        self._material = append_material("MN Flat Outline")


class TransparentOutline(MaterialTreeInterface):
    def __init__(self):
        self._material = append_material("MN Transparent Outline")


class Squishy(MaterialTreeInterface):
    def __init__(self):
        self._material = append_material("MN Squishy")

    subsurface_scale = socket("Principled BSDF", "Subsurface Scale", float)


class Default(MaterialTreeInterface):
    def __init__(self):
        self._material = append_material("MN Default")

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
