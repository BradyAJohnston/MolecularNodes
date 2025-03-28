import bpy
import numpy as np
from bpy.types import Material, ShaderNodeTree
from databpy.material import append_from_blend

from ..nodes import TreeInterface
from ..interface import option, socket
from ...assets import MN_DATA_FILE


MATERIAL_NAMES = [
    "MN Default",
    "MN Flat Outline",
    "MN Squishy",
    "MN Transparent Outline",
    "MN Ambient Occlusion",
]


def append_material(name: str) -> bpy.types.Material:
    "Append a material from the MN_DATA_FILE."
    return append_from_blend(name, str(MN_DATA_FILE))


def add_all_materials() -> None:
    "Append all pre-defined materials from the MN_DATA_FILE."
    for name in MATERIAL_NAMES:
        append_material(name)


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

    def _expose_all_inputs(self):
        for node in self.node_tree.nodes:
            for input in node.inputs:
                if hasattr(input, "default_value"):
                    auto_socket(input)


def dynamic_tree_interface(tree: bpy.types.ShaderNodeTree) -> MaterialTreeInterface:
    class_name = f'DynamicMaterialInterface_{tree.name.replace(" ", "_")}'

    DynaicMaterialInterface = type(
        class_name,
        (MaterialTreeInterface,),
        {"_material": None},
    )

    interface = DynaicMaterialInterface._from_bpy_material(tree)
    interface._expose_all_inputs()
    return interface
