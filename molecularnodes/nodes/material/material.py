import bpy
from bpy.types import Material, ShaderNodeTree
from databpy.material import append_from_blend

from ..interface import option, socket, TreeInterface
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
    def __init__(self, material: Material):
        if isinstance(material, str):
            material = bpy.data.materials[material]
        elif isinstance(material, Material):
            material = material
        else:
            raise ValueError("Material must be a string or a Material object")
        self.material = material
        self.material: Material = material

    @property
    def node_tree(self) -> ShaderNodeTree:
        if self.material.node_tree is None:
            raise ValueError("Material has no node tree")
        return self.material.node_tree

    @property
    def nodes(self):
        return self.node_tree.nodes

    @property
    def links(self):
        return self.node_tree.links

    def _expose_all_inputs(self):
        for node in self.node_tree.nodes:
            if "Material Output" in node.name:
                continue
            for input in node.inputs:
                prop_name = input.name.lower().replace(" ", "_")
                if hasattr(input, "default_value"):
                    setattr(self.__class__, prop_name, socket(input))


def assign_material(node, new_material: str | bpy.types.Material = "default") -> None:
    add_all_materials()
    material_socket = node.inputs.get("Material")
    if material_socket is None:
        return None

    if isinstance(new_material, bpy.types.Material):
        material_socket.default_value = new_material
    elif new_material == "default":
        material_socket.default_value = append_material("MN Default")
    else:
        try:
            material_socket.default_value = append_material(new_material)
        except Exception as e:
            print(f"Unable to use material {new_material}, error: {e}")


def dynamic_material_interface(material: bpy.types.Material) -> MaterialTreeInterface:
    class_name = (
        f'DynamicMaterialInterface_{material.name.replace(" ", "_").replace(".", "_")}'
    )

    DynaicMaterialInterface = type(
        class_name,
        (MaterialTreeInterface,),
        {},
    )

    interface = DynaicMaterialInterface(material)
    interface._expose_all_inputs()
    return interface
