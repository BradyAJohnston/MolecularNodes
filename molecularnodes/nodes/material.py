import bpy
from bpy.types import Material, ShaderNodeTree
from databpy.material import append_from_blend
from ..assets import MN_DATA_FILE
from .interface import TreeInterface, check_linked, remove_linked, socket

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


# class to interact with a bpy.types.Material node tree and change some of the default
# values of the nodes inside of it
class MaterialTreeInterface(TreeInterface):
    def __init__(self, material: Material):
        super().__init__()
        self._allowable_properties.add("material")
        if isinstance(material, str):
            material = bpy.data.materials[material]
        elif isinstance(material, Material):
            material = material
        else:
            raise ValueError("Material must be a string or a Material object")

        self.material: Material = material
        self.tree: ShaderNodeTree = self.material.node_tree  # type: ignore

    def _expose_all_inputs(self):
        for node in self.tree.nodes:
            if "Material Output" in node.name:
                continue
            for input in node.inputs:
                if input.is_unavailable:
                    continue
                input_name = input.name.lower().replace(" ", "_")
                node_name = (
                    (node.name if node.label == "" else node.label)
                    .lower()
                    .replace(" ", "_")
                )
                prop_name = f"{node_name}_{input_name}"

                if hasattr(input, "default_value"):
                    self._register_property(prop_name)
                    setattr(self.__class__, prop_name, socket(input))


def set_socket_material(
    socket: bpy.types.NodeSocketMaterial,
    mat: MaterialTreeInterface
    | bpy.types.Material
    | bpy.types.NodeSocketMaterial
    | str
    | None,
) -> None:
    remove_linked(socket)
    if mat is None:
        socket.default_value = None
    elif isinstance(mat, MaterialTreeInterface):
        socket.default_value = mat.material
    elif isinstance(mat, bpy.types.Material):
        socket.default_value = mat
    elif isinstance(mat, bpy.types.NodeSocketMaterial):
        socket.node.id_data.links.new(mat, socket)  # type: ignore
    elif isinstance(mat, str):
        mat = bpy.data.materials[mat]
        socket.default_value = mat
    else:
        raise TypeError("Invalid type for setting of a material: " + str(type(mat)))


def assign_material(
    node: bpy.types.GeometryNodeGroup,
    new_material: MaterialTreeInterface
    | bpy.types.Material
    | bpy.types.NodeSocketMaterial
    | str
    | None = "default",
) -> None:
    add_all_materials()
    if isinstance(new_material, str):
        if new_material not in bpy.data.materials:
            try:
                append_material(new_material)
            except Exception:
                try:
                    new_material = "MN " + new_material.title().strip()
                    append_material(new_material)
                except Exception:
                    raise ValueError(
                        f"Material {new_material} not found in this file of the included MN preset file."
                    )
    try:
        set_socket_material(
            socket=node.inputs["Material"],
            mat=new_material,
        )
    except KeyError:
        pass


def dynamic_material_interface(material: bpy.types.Material) -> MaterialTreeInterface:
    class_name = (
        f"DynamicMaterialInterface_{material.name.replace(' ', '_').replace('.', '_')}"
    )

    DynaicMaterialInterface = type(
        class_name,
        (MaterialTreeInterface,),
        {},
    )

    interface = DynaicMaterialInterface(material)
    interface._expose_all_inputs()
    return interface


def getset_material(socket: bpy.types.NodeSocketMaterial):
    def getter(self) -> MaterialTreeInterface | None:
        check_linked(socket)
        mat = getattr(socket, "default_value")
        if mat is None:
            return None
        else:
            interface = dynamic_material_interface(mat)
            return interface

    def setter(
        self,
        mat: MaterialTreeInterface
        | bpy.types.Material
        | bpy.types.NodeSocketMaterial
        | str
        | None,
    ) -> None:
        set_socket_material(socket, mat)

    return property(getter, setter)


class MaterialConstructor(MaterialTreeInterface):
    def __init__(self, material_name: str, **kwargs):
        super().__init__(append_material(material_name))
        self._expose_all_inputs()
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                raise ValueError(f"Material does not have property {key}")


class AmbientOcclusion(MaterialConstructor):
    def __init__(self, **kwargs):
        super().__init__("MN Ambient Occlusion", **kwargs)


class Default(MaterialConstructor):
    def __init__(self, **kwargs):
        super().__init__("MN Default", **kwargs)


class FlatOutline(MaterialConstructor):
    def __init__(self, **kwargs):
        super().__init__("MN Flat Outline", **kwargs)


class Squishy(MaterialConstructor):
    def __init__(self, **kwargs):
        super().__init__("MN Squishy", **kwargs)


class TransparentOutline(MaterialConstructor):
    def __init__(self, **kwargs):
        super().__init__("MN Transparent Outline", **kwargs)
