from typing import Tuple
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
        elif isinstance(
            material, bpy.types.Material
        ):  # note this is equivalent class to the above.
            material = material
        else:
            raise ValueError(
                f"Material must be a string or a Material object. Found {type(material)}"
            )

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
        return "Material input not found on node."


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


class Material:
    """
    Note: if this is a nice route we could add fns that will create MN_materials here and remove them from the blend file.

    Materials:
    - AmbientOcclusion
    - Default  ( tried to replicate the defaults below )
    - FlatOutlineq
    - Squishy
    - TransparentOutline

    """

    @staticmethod
    def glass() -> bpy.types.Material:
        return create_material(
            name="MN Glass",  # Add a descriptive name
            base_color=(0.8, 0.9, 1.0, 0.2),
            transmission_weight=0.95,
            roughness=0.0,
            ior=1.45,
            specular_ior_level=0.6,  # Note: Renamed to "Specular IOR" in create_material
        )

    @staticmethod
    def gold() -> bpy.types.Material:
        return create_material(
            name="MN Gold",
            base_color=(1.0, 0.8, 0.2, 1.0),
            metallic=1.0,
            roughness=0.3,
            anisotropic=0.8,
            anisotropic_rotation=0.5,
        )

    @staticmethod
    def green_glow() -> bpy.types.Material:
        return create_material(
            name="MN GreenGlow",
            base_color=(0.0, 1.0, 0.0, 1.0),
            emission_strength=5.0,
            emission_color=(0.0, 1.0, 0.0, 1.0),
            metallic=0.0,
            roughness=0.2,
        )

    @staticmethod
    def holo() -> bpy.types.Material:
        return create_material(
            name="MN Holo",
            base_color=(0.2, 0.6, 1.0, 0.3),
            emission_strength=2.0,
            emission_color=(0.2, 0.6, 1.0, 1.0),
            transmission_weight=0.8,
            thin_film_thickness=1000,
        )

    @staticmethod
    def iridescent() -> bpy.types.Material:
        return create_material(
            name="MN Iridescent",
            base_color=(0.8, 0.8, 0.8, 0.3),
            metallic=0.8,
            transmission_weight=0.5,
            thin_film_thickness=1200,
            thin_film_ior=1.5,
            coat_weight=1.0,
        )

    @staticmethod
    def metallic() -> bpy.types.Material:
        return create_material(
            name="MN Metallic",
            base_color=(0.7, 0.7, 0.7, 1.0),
            metallic=1.0,
            roughness=0.1,
            specular_ior_level=0.9,
            coat_weight=1.0,
        )

    @staticmethod
    def neon() -> bpy.types.Material:
        return create_material(
            name="MN Neon",
            base_color=(0.0, 1.0, 0.8, 1.0),
            emission_strength=5.0,
            emission_color=(0.0, 1.0, 0.8, 1.0),
            metallic=0.8,
            roughness=0.1,
        )

    @staticmethod
    def new() -> bpy.types.Material:
        return create_material(name="MN Default 02")

    @staticmethod
    def pearl() -> bpy.types.Material:
        return create_material(
            name="MN Pearl",
            base_color=(0.9, 0.9, 0.9, 1.0),
            metallic=0.7,
            roughness=0.15,
            coat_weight=1.0,
            coat_ior=2.0,
            thin_film_thickness=500,
        )

    @staticmethod
    def subsurface() -> bpy.types.Material:
        return create_material(
            name="MN Subsurface",
            base_color=(1.0, 0.4, 0.4, 1.0),
            subsurface_weight=1.0,
            subsurface_radius=(1.0, 0.2, 0.1),
            subsurface_scale=0.5,
            emission_strength=0.3,
        )

    @staticmethod
    def toon() -> bpy.types.Material:
        return create_material(
            name="MN Toon",
            base_color=(0.2, 0.6, 1.0, 1.0),
            metallic=0,
            roughness=1.0,
            specular_ior_level=0.0,
            diffuse_roughness=1.0,
            coat_weight=0.2,
        )

    @staticmethod
    def velvet() -> bpy.types.Material:
        return create_material(
            name="MN Velvet",
            base_color=(0.5, 0.0, 0.2, 1.0),
            sheen_weight=1.0,
            sheen_roughness=0.3,
            sheen_tint=(1.0, 0.8, 0.9, 1.0),
            roughness=0.8,
        )

    @staticmethod
    def waxy() -> bpy.types.Material:
        return create_material(
            name="MN Waxy",
            base_color=(0.9, 0.87, 0.82, 1.0),
            subsurface_weight=0.3,
            subsurface_radius=(0.5, 0.4, 0.3),
            subsurface_scale=0.1,
            roughness=0.4,
            specular_ior_level=0.2,
            metallic=0.0,
            coat_weight=0.1,
            coat_roughness=0.3,
            sheen_weight=0.1,
            sheen_roughness=0.3,
        )


def create_material(
    name: str | None = None,
    base_color: Tuple[float, float, float, float] = (0.8, 0.8, 0.8, 0.05),
    metallic: float = 0.0,
    roughness: float = 0.2,
    ior: float = 1.45,
    alpha: float = 1.0,
    normal: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    weight: float = 0.0,
    diffuse_roughness: float = 0.0,
    subsurface_weight: float = 0.0,
    subsurface_radius: Tuple[float, float, float] = (1.0, 0.2, 0.1),
    subsurface_scale: float = 0.05,
    subsurface_ior: float = 1.4,
    subsurface_anisotropy: float = 0.0,
    specular_ior_level: float = 0.5,
    specular_tint: Tuple[float, float, float, float] = (1.0, 1.0, 1.0, 1.0),
    anisotropic: float = 0.0,
    anisotropic_rotation: float = 0.0,
    tangent: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    transmission_weight: float = 0.0,
    coat_weight: float = 0.0,
    coat_roughness: float = 0.03,
    coat_ior: float = 1.5,
    coat_tint: Tuple[float, float, float, float] = (1.0, 1.0, 1.0, 1.0),
    coat_normal: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    sheen_weight: float = 0.0,
    sheen_roughness: float = 0.5,
    sheen_tint: Tuple[float, float, float, float] = (0.5, 0.5, 0.5, 1.0),
    emission_color: Tuple[float, float, float, float] = (0.0, 0.0, 0.0, 1.0),
    emission_strength: float = 0.0,
    thin_film_thickness: float = 0.0,
    thin_film_ior: float = 1.3,
) -> bpy.types.Material:
    """Create a Blender material from provided shader properties."""

    key_map = {
        "Base Color": base_color,
        "Metallic": metallic,
        "Roughness": roughness,
        "IOR": ior,
        "Alpha": alpha,
        "Normal": normal,
        "Weight": weight,
        "Diffuse Roughness": diffuse_roughness,
        "Subsurface Weight": subsurface_weight,
        "Subsurface Radius": subsurface_radius,
        "Subsurface Scale": subsurface_scale,
        "Subsurface IOR": subsurface_ior,
        "Subsurface Anisotropy": subsurface_anisotropy,
        "Specular IOR Level": specular_ior_level,
        "Specular Tint": specular_tint,
        "Anisotropic": anisotropic,
        "Anisotropic Rotation": anisotropic_rotation,
        "Tangent": tangent,
        "Transmission Weight": transmission_weight,
        "Coat Weight": coat_weight,
        "Coat Roughness": coat_roughness,
        "Coat IOR": coat_ior,
        "Coat Tint": coat_tint,
        "Coat Normal": coat_normal,
        "Sheen Weight": sheen_weight,
        "Sheen Roughness": sheen_roughness,
        "Sheen Tint": sheen_tint,
        "Emission Color": emission_color,
        "Emission Strength": emission_strength,
        "Thin Film Thickness": thin_film_thickness,
        "Thin Film IOR": thin_film_ior,
    }

    if name is None:
        name = f"material_func_{id(key_map)}"

    if name in bpy.data.materials:
        return bpy.data.materials[name]

    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes.get("Principled BSDF")

    for input_socket in bsdf.inputs:
        if input_socket.name in key_map:
            input_socket.default_value = key_map[input_socket.name]

    return mat
