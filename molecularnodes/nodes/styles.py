"""
Style Classes

MolecularNodes uses Geometry Nodes to define molecule representations within blender. However,
it is desireable to make these nodes accessible for scripting. The classes defined in this file
provide a way to script against the GeometryNode nodetrees. When invoked on a style these classes
will simply override the base styles.

Because the nodetrees in the MN_data_file_{version}.blend are the source of truth, we add
function to parity between the nodetrees and the class representations.

"""

from typing import List, Tuple, Union, Any, Dict
from bpy.types import GeometryNodeGroup
# from .nodes import styles_mapping
# import databpy import db
# from databpy.nodes import (append_from_blend)

__all__ = [
    "StyleBallandStick",
    "StyleCartoon",
    "StyleRibbon",
    "StyleSpheres",
    "StyleSticks",
    "StyleSurface",
    "StyleBase",
]


PortDataEntry = Dict[str, Any]
PortDataList = List[PortDataEntry]


class StyleBase:
    style: str
    portdata: PortDataList = []

    def update_style_node(self, node_style: GeometryNodeGroup):
        for input in node_style.inputs:
            if input.type != "GEOMETRY":
                for arg in self.portdata:
                    name = arg["name"]
                    blendername = arg.get("blendername")
                    if input.name == blendername:
                        input.default_value = getattr(self, name)


class StyleBallandStick(StyleBase):
    style = "ball_and_stick"
    portdata: PortDataList = [
        {"name": "quality", "blendername": "Quality"},
        {"name": "geometry", "blendername": "Sphere Geometry"},
        {"name": "sphere_radii", "blendername": "Sphere Radii"},
        {"name": "bond_split", "blendername": "Bond Split"},
        {"name": "bond_find", "blendername": "Bond Find"},
        {"name": "bond_radius", "blendername": "Bond Radius"},
        {"name": "color_blur", "blendername": "Color Blur"},
        {"name": "shade_smooth", "blendername": "Shade Smooth"},
    ]

    def __init__(
        self,
        quality: int = 2,
        geometry: str = "Instance",  # enum
        sphere_radii: float = 0.3,
        bond_split: bool = False,
        bond_find: bool = False,
        bond_radius: float = 0.3,
        color_blur: bool = False,
        shade_smooth: bool = True,
    ):
        self.quality = quality
        self.geometry = geometry
        self.sphere_radii = sphere_radii
        self.bond_split = bond_split
        self.bond_find = bond_find
        self.bond_radius = bond_radius
        self.color_blur = color_blur
        self.shade_smooth = shade_smooth


class StyleCartoon(StyleBase):
    style = "cartoon"
    portdata: PortDataList = [
        {"name": "quality", "blendername": "Quality"},
        {"name": "dssp", "blendername": "Peptide DSSP"},
        {"name": "cylinders", "blendername": "Peptide Cylinders"},
        {"name": "arrows", "blendername": "Peptide Arrows"},
        {"name": "rounded", "blendername": "Peptide Rounded"},
        {"name": "thickness", "blendername": "Peptide Thickness"},
        {"name": "width", "blendername": "Peptide Width"},
        {"name": "loop_radius", "blendername": "Peptide Loop Radius"},
        {"name": "smoothing", "blendername": "Peptide Smoothing"},
        {"name": "backbone_shape", "blendername": "Backbone Shape"},
        {"name": "backbone_radius", "blendername": "Backbone Radius"},
        {"name": "base_shape", "blendername": "Base Shape"},
        {"name": "base_realize", "blendername": "Base Realize"},
        {"name": "color_blur", "blendername": "Color Blur"},
        {"name": "shade_smooth", "blendername": "Shade Smooth"},
    ]

    def __init__(
        self,
        quality: int = 2,
        dssp: bool = False,
        cylinders: bool = False,
        arrows: bool = True,
        rounded: bool = False,
        thickness: float = 0.6,
        width: float = 2.2,
        loop_radius: float = 0.3,
        smoothing: float = 0.5,
        backbone_shape: str = "Cylinder",  # enum?
        backbone_radius: float = 2.0,
        base_shape: str = "Rectangle",  # enum?
        base_realize: bool = False,
        color_blur: bool = True,
        shade_smooth: bool = True,
    ):
        self.quality = quality
        self.dssp = dssp
        self.cylinders = cylinders
        self.arrows = arrows
        self.rounded = rounded
        self.thickness = thickness
        self.width = width
        self.loop_radius = loop_radius
        self.smoothing = smoothing
        self.backbone_shape = backbone_shape
        self.backbone_radius = backbone_radius
        self.base_shape = base_shape
        self.base_realize = base_realize
        self.color_blur = color_blur
        self.shade_smooth = shade_smooth


class StyleRibbon(StyleBase):
    style = "ribbon"
    portdata: PortDataList = [
        {"name": "quality", "blendername": "Quality"},
        {"name": "radius", "blendername": "Radius"},
        {"name": "smoothing", "blendername": "Smoothing"},
        {"name": "color_blur", "blendername": "Color Blur"},
        {"name": "shade_smooth", "blendername": "Shade Smooth"},
        {"name": "backbone_smoothing", "blendername": "Backbone Smoothing"},
        {"name": "backbone_threshold", "blendername": "Backbone Threshold"},
        {"name": "backbone_radius", "blendername": "Backbone Radius"},
        {"name": "backbone_shape", "blendername": "Backbone Shape"},
        {"name": "base_scale", "blendername": "Base Scale"},
        {"name": "base_resolution", "blendername": "Base Resolution"},
        {"name": "base_realize", "blendername": "Base Realize"},
        {"name": "uv_map", "blendername": "UV Map"},
        {"name": "u_component", "blendername": "U Component"},
        {"name": "u_component_factor", "blendername": "U Component Factor"},
    ]

    def __init__(
        self,
        quality: int = 3,
        radius: float = 1.6,
        smoothing: float = 0.6,
        color_blur: bool = False,
        shade_smooth: bool = False,
        backbone_smoothing: float = 0.5,
        backbone_threshold: float = 4.5,
        base_scale: Tuple[float, float, float] = (2.5, 0.5, 7.0),
        backbone_radius: float = 1.6,
        # backbone_shape: str = "Cylinder",
        base_resolution: int = 4,
        base_realize: bool = False,
        uv_map: bool = False,
        u_component=None,
        u_component_factor=None,
    ):
        self.quality = quality
        self.radius = radius
        self.smoothing = smoothing
        self.color_blur = color_blur
        self.shade_smooth = shade_smooth
        self.backbone_smoothing = backbone_smoothing
        self.backbone_threshold = backbone_threshold
        self.backbone_radius = backbone_radius
        # self.backbone_shape = backbone_shape  # enum?
        self.base_scale = base_scale
        self.base_resolution = base_resolution
        self.base_realize = base_realize
        self.uv_map = uv_map
        self.u_component = u_component
        self.u_component_factor = u_component_factor


class StyleSpheres(StyleBase):
    style = "spheres"
    portdata: PortDataList = [
        {"name": "geometry", "blendername": "Sphere Geometry"},
        {"name": "radii", "blendername": "Sphere Radii"},
        {"name": "sphere_subdivisions", "blendername": "Sphere Subdivisions"},
        {"name": "shade_smooth", "blendername": "Shade Smooth"},
    ]

    def __init__(
        self,
        geometry: str = "Point",  # enum: "Point" (Point Cloud), "Instances" (Instances of a mesh Icosphere), or "Mesh" (realised Mesh)
        radii: float = 0.8,
        sphere_subdivisions: int = 2,
        shade_smooth: bool = True,
    ):
        self.geometry = geometry
        self.radii = radii
        self.sphere_subdivisions = sphere_subdivisions
        self.shade_smooth = shade_smooth


class StyleSticks(StyleBase):
    style = "sticks"
    portdata: PortDataList = [
        {"name": "quality", "blendername": "Quality"},
        {"name": "radius", "blendername": "Radius"},
        {"name": "color_blur", "blendername": "Color Blur"},
        {"name": "shade_smooth", "blendername": "Shade Smooth"},
    ]

    def __init__(
        self,
        quality: int = 3,
        radius: float = 0.20000000298023224,
        color_blur: bool = False,
        shade_smooth: bool = True,
    ):
        self.quality = quality
        self.radius = radius
        self.color_blur = color_blur
        self.shade_smooth = shade_smooth


class StyleSurface(StyleBase):
    style = "surface"
    portdata: PortDataList = [
        {"name": "quality", "blendername": "Quality"},
        {"name": "scale_radii", "blendername": "Scale Radii"},
        {"name": "probe_size", "blendername": "Probe Size"},
        {"name": "relaxation_steps", "blendername": "Relaxation Steps"},
        {"name": "separate", "blendername": "Separate By"},
        {"name": "group_id", "blendername": "Group ID"},
        {"name": "color_source", "blendername": "Color Source"},
        {"name": "blur", "blendername": "Color Blur"},
        {"name": "shade_smooth", "blendername": "Shade Smooth"},
    ]

    def __init__(
        self,
        quality: int = 3,
        scale_radii: float = 1.5,
        probe_size: float = 1.0,
        relaxation_steps: int = 10,
        separate: str = "chain_id",  # enum?
        group_id: int = 0,
        color_source: str = "Alpha Carbon",  # enum?
        blur: int = 2,
        shade_smooth: bool = True,
    ):
        self.quality = quality
        self.scale_radii = scale_radii
        self.probe_size = probe_size
        self.relaxation_steps = relaxation_steps
        self.separate = separate
        self.group_id = group_id
        self.color_source = color_source
        self.blur = blur
        self.shade_smooth = shade_smooth
