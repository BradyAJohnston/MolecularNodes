import bpy
from databpy.material import append_from_blend
from ..assets import MN_DATA_FILE
from .interface import remove_linked

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


def add_all_materials() -> dict[str, bpy.types.Material]:
    "Append all pre-defined materials from the MN_DATA_FILE."
    return {name: append_material(name) for name in MATERIAL_NAMES}


def set_socket_material(
    socket: bpy.types.NodeSocketMaterial,
    mat: bpy.types.Material | bpy.types.NodeSocketMaterial | str | None,
) -> None:
    remove_linked(socket)
    if mat is None:
        socket.default_value = None
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
    new_material: bpy.types.Material
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
        return "Material input not found on node."
