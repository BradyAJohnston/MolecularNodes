import bpy
import numpy as np
from bpy.types import Material, ShaderNodeTree
from typing import Literal

from .nodes import TreeInterface
from .utils import append_material, option, socket


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

    @classmethod
    def _from_bpy_material(cls, material: str | Material):
        mat = cls()
        if isinstance(material, str):
            material = bpy.data.materials[material]
        elif isinstance(material, Material):
            material = material
        else:
            raise ValueError("Material must be a string or a Material object")
        mat._material = material
        return mat


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


def generic_material_interface(material: Material) -> MaterialTreeInterface:
    # Create a new class that inherits from MaterialTreeInterface
    class GenericInterface(MaterialTreeInterface):
        def __init__(self):
            self._material = material

    # For each node in the material's node tree
    for node in material.node_tree.nodes:  # type: ignore
        for input in node.inputs:
            # Only process unlinked inputs
            if input.is_linked:
                continue
            try:
                value = input.default_value  # type: ignore
            except AttributeError:
                continue

            # Determine the type based on the value
            if isinstance(value, (int, float, bool, str)):
                value_type = type(value)
            else:
                value_type = np.ndarray
                value = np.array(value)

            # Create property name from node and input names
            prop_name = f"{node.name}_{input.name}".lower().replace(" ", "_")

            # Add the socket as a property to the class
            # shader sockets aren't to be accessed directly so just ignore them
            try:
                setattr(
                    GenericInterface,
                    prop_name,
                    socket(node.name, input.name, value_type),  # type: ignore
                )
            except AttributeError:
                pass

    return GenericInterface()


class AmbientOcclusion(MaterialTreeInterface):
    def __init__(self):
        self._material = append_material("MN Ambient Occlusion")

    amount = socket("Mix", "Factor", float)
    power = socket("Math", 1, float)
    distance = socket("Ambient Occlusion", "Distance", float)
    samples = option("Ambient Occlusion", "samples", int)
    inside = option("Ambient Occlusion", "inside", bool)
    only_local = option("Ambient Occlusion", "Only Local", bool)


class FlatOutline(MaterialTreeInterface):
    def __init__(self):
        self._material = append_material("MN Flat Outline")


class TransparentOutline(MaterialTreeInterface):
    def __init__(self):
        self._material = append_material("MN Transparent Outline")


class Default(MaterialTreeInterface):
    def __init__(self):
        self._material = append_material("MN Default")

    metallic = socket("Principled BSDF", "Metallic", float)
    roughness = socket("Principled BSDF", "Roughness", float)
    ior = socket("Principled BSDF", "IOR", float)
    alpha = socket("Principled BSDF", "Alpha", float)

    diffuse_roughness = socket("Principled BSDF", "Diffuse Roughness", float)

    subsurface_method = option(
        "Principled BSDF",
        "Subsurface Method",
        str,
    )
    subsurface_weight = socket("Principled BSDF", "Subsurface Weight", float)
    subsurface_radius = socket(
        "Principled BSDF", "Subsurface Radius", tuple[float, float, float]
    )
    subsurface_scale = socket("Principled BSDF", "Subsurface Scale", float)


class Squishy(Default):
    def __init__(self):
        self._material = append_material("MN Squishy")
