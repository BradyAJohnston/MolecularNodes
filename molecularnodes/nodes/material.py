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
        return f"Material input not found on node."


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


import bpy
from dataclasses import dataclass, replace, field, fields
from typing import List, Tuple, Union


__all__ = ['Material', 'MaterialCreator']


@dataclass(frozen=True)
class Material():
    """Provides all material properties available to a Blender BSDF Material.

    For full details on Blender material properties, see the Blender documentation:
    https://docs.blender.org/manual/en/latest/render/shader_nodes/shader/principled.html

    Material instances store all properties for Blender's Principled BSDF shader node.
    Typically created using MaterialCreator's factory methods or the Material.update_properties() method.

    Examples:
        >>> # Create a basic material
        >>> mat = Material()
        >>> # Update properties of an existing material
        >>> transparent_mat = mat.update_properties(alpha=0.5)
    """
    base_color: Tuple[float, float, float, float] = field(default=(0.8, 0.8, 0.8, 0.05), metadata={"key": "Base Color"})
    metallic: float = field(default=0.0, metadata={"key": "Metallic"})
    roughness: float = field(default=0.2, metadata={"key": "Roughness"})
    ior: float = field(default=1.45, metadata={"key": "IOR"})
    alpha: float = field(default=1.0, metadata={"key": "Alpha"})
    normal: Tuple[float, float, float] = field(default=(0.0, 0.0, 0.0), metadata={"key": "Normal"})
    weight: float = field(default=0.0, metadata={"key": "Weight"})
    diffuse_roughness: float = field(default=0.0, metadata={"key": "Diffuse Roughness"})
    subsurface_weight: float = field(default=0.0, metadata={"key": "Subsurface Weight"})
    subsurface_radius: Tuple[float, float, float] = field(default=(1.0, 0.2, 0.1), metadata={"key": "Subsurface Radius"})
    subsurface_scale: float = field(default=0.05, metadata={"key": "Subsurface Scale"})
    subsurface_ior: float = field(default=1.4, metadata={"key": "Subsurface IOR"})
    subsurface_anisotropy: float = field(default=0.0, metadata={"key": "Subsurface Anisotropy"})
    specular_ior_level: float = field(default=0.5, metadata={"key": "Specular IOR Level"})
    specular_tint: Tuple[float, float, float, float] = field(default=(1.0, 1.0, 1.0, 1.0), metadata={"key": "Specular Tint"})
    anisotropic: float = field(default=0.0, metadata={"key": "Anisotropic"})
    anisotropic_rotation: float = field(default=0.0, metadata={"key": "Anisotropic Rotation"})
    tangent: Tuple[float, float, float] = field(default=(0.0, 0.0, 0.0), metadata={"key": "Tangent"})
    transmission_weight: float = field(default=0.0, metadata={"key": "Transmission Weight"})
    coat_weight: float = field(default=0.0, metadata={"key": "Coat Weight"})
    coat_roughness: float = field(default=0.03, metadata={"key": "Coat Roughness"})
    coat_ior: float = field(default=1.5, metadata={"key": "Coat IOR"})
    coat_tint: Tuple[float, float, float, float] = field(default=(1.0, 1.0, 1.0, 1.0), metadata={"key": "Coat Tint"})
    coat_normal: Tuple[float, float, float] = field(default=(0.0, 0.0, 0.0), metadata={"key": "Coat Normal"})
    sheen_weight: float = field(default=0.0, metadata={"key": "Sheen Weight"})
    sheen_roughness: float = field(default=0.5, metadata={"key": "Sheen Roughness"})
    sheen_tint: Tuple[float, float, float, float] = field(default=(0.5, 0.5, 0.5, 1.0), metadata={"key": "Sheen Tint"})
    emission_color: Tuple[float, float, float, float] = field(default=(0.0, 0.0, 0.0, 1.0), metadata={"key": "Emission Color"})
    emission_strength: float = field(default=0.0, metadata={"key": "Emission Strength"})
    thin_film_thickness: float = field(default=0.0, metadata={"key": "Thin Film Thickness"})
    thin_film_ior: float = field(default=1.3, metadata={"key": "Thin Film IOR"})

    def get_by_key(self, original_key: str):
        for f in fields(self):
            if f.metadata.get('key') == original_key:
                return getattr(self, f.name)
        return None

    def update_properties(self, **changes):
        return replace(self, **changes)

    def blenderize(self, name=None) -> bpy.types.Material:
        """Convert this Material to a Blender material.

         Args:
             name (str, optional): Name for the Blender material. Defaults to a generated unique name.

         Returns:
             bpy.types.Material: The created or retrieved Blender material

         Examples:
             >>> glass_mat = MaterialCreator.glass()
             >>> blender_material = glass_mat.blenderize("MyGlassMaterial")
        """
        if name is None:
            name = f"material_{id(self)}"

        if name in bpy.data.materials:
            return bpy.data.materials[name]
        else:
            mat = bpy.data.materials.new(name)
            mat.use_nodes = True
            bsdf = mat.node_tree.nodes.get("Principled BSDF")
            for input in bsdf.inputs:
                if input.type != "GEOMETRY":
                    if value := self.get_by_key(input.name):
                        input.default_value = value
            return mat


class MaterialCreator():
    """Factory class that provides pre-configured material presets.

    This class contains static methods to create various material presets for use
    in Blender. Each method returns a Material instance configured with specific
    properties to achieve common material effects.

    Examples:
        >>> # Create a glass material
        >>> glass_material = MaterialCreator.glass()
        >>> # Create a metal material and apply it to an object in Blender
        >>> metal_material = MaterialCreator.metallic()
        >>> blender_mat = metal_material.blenderize("MetalMaterial")
        >>> my_object.data.materials.append(blender_mat)
    """
    def __init__(self):
        pass

    def glass() -> Material:
        return replace(Material(),
            base_color=(0.8, 0.9, 1.0, 0.2),
            transmission_weight=0.95,
            roughness=0.0,
            ior=1.45,
            specular_ior_level=0.6
        )

    def gold() -> Material:
        return replace(Material(),
            base_color=(1.0, 0.8, 0.2, 1.0),
            metallic=1.0,
            roughness=0.3,
            anisotropic=0.8,
            anisotropic_rotation=0.5
        )

    def green_glow() -> Material:
        return replace(Material(),
            base_color=(0.0, 1.0, 0.0, 1.0),
            emission_strength=5.0,
            emission_color=(0.0, 1.0, 0.0, 1.0),
            metallic=0.0,
            roughness=0.2
            )

    def holo() -> Material:
        return replace(Material(),
            base_color=(0.2, 0.6, 1.0, 0.3),
            emission_strength=2.0,
            emission_color=(0.2, 0.6, 1.0, 1.0),
            transmission_weight=0.8,
            thin_film_thickness=1000
        )

    def iridescent() -> Material:
        return replace(Material(),
            base_color=(0.8, 0.8, 0.8, 0.3),
            metallic=0.8,
            transmission_weight=0.5,
            thin_film_thickness=1200,
            thin_film_ior=1.5,
            coat_weight=1.0
        )

    def metallic() -> Material:
        return replace(Material(),
            base_color=(0.7, 0.7, 0.7, 1.0),
            metallic=1.0,
            roughness=0.1,
            specular_ior_level=0.9,
            coat_weight=1.0
        )

    def neon() -> Material:
        return replace(Material(),
            base_color=(0.0, 1.0, 0.8, 1.0),
            emission_strength=5.0,
            emission_color=(0.0, 1.0, 0.8, 1.0),
            metallic=0.8,
            roughness=0.1)

    def new() -> Material:
        return Material()

    def pearl() -> Material:
        return replace(Material(),
            base_color=(0.9, 0.9, 0.9, 1.0),
            metallic=0.7,
            roughness=0.15,
            coat_weight=1.0,
            coat_ior=2.0,
            thin_film_thickness=500
        )

    def subsurface() -> Material:
        return replace(Material(),
            base_color=(1.0, 0.4, 0.4, 1.0),
            subsurface_weight=1.0,
            subsurface_radius=(1.0, 0.2, 0.1),
            subsurface_scale=0.5,
            emission_strength=0.3
        )

    def toon() -> Material:
        return replace(Material(),
            base_color=(0.2, 0.6, 1.0, 1.0),
            metallic=0,
            roughness=1.0,
            specular_ior_level=0.0,
            diffuse_roughness=1.0,
            coat_weight=0.2
        )

    def velvet() -> Material:
        return replace(Material(),
            base_color=(0.5, 0.0, 0.2, 1.0),
            sheen_weight=1.0,
            sheen_roughness=0.3,
            sheen_tint=(1.0, 0.8, 0.9, 1.0),
            roughness=0.8
        )

    def waxy() -> Material:
        return replace(Material(),
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
            sheen_roughness=0.3
        )
