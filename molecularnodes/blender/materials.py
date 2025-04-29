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
