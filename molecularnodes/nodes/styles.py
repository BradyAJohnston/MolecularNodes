from dataclasses import dataclass, replace, field, fields
from typing import List, Tuple, Union, Any, Dict
from bpy.types import GeometryNodeGroup

__all__ = [
    "StyleBallandStick",
    "StyleCartoon",
    "StyleRibbon",
    "StyleSpheres",
    "StyleSticks",
    "StyleSurface",
    "StyleBase",
]


# Define the type for a single port data entry for clarity
PortDataEntry = Dict[str, Any]
# Define the type for the list of port data entries
PortDataList = List[PortDataEntry]


class StyleBase:
    portdata: PortDataList = []
    def update_style_node(self, node_style: GeometryNodeGroup):
        for input in node_style.inputs:
            if input.type != "GEOMETRY":
                print(input.name, input.default_value,)
                for arg in self.portdata:
                    # print(arg)
                    name = arg['name']
                    blendername = arg.get('blendername', name)  # Use name if blendername not
                    if input.name == blendername:
                        #print(f"setting! Val : {getattr(self, name)}")
                        #print(f"Initial Val : {input.default_value}")
                        input.default_value = getattr(self, name)
                        #print(f"AfterSetting Val : {input.default_value}")


class StyleBallandStick(StyleBase):
    # fmt: off
    portdata: PortDataList = [
        { "name": "quality", "blendername": "Quality", "type": int, "default": 2 },
        { "name": "sphere_radii", "blendername": "Sphere Radii", "type": float, "default": 0.30000001192092896 },
        { "name": "bond_split", "blendername": "Bond Split", "type": bool, "default": False },
        { "name": "bond_find", "blendername": "Bond Find", "type": bool, "default": False },
        { "name": "bond_radius", "blendername": "Bond Radius", "type": float, "default": 0.30000001192092896 },
        { "name": "color_blur", "blendername": "Color Blur", "type": bool, "default": False },
        { "name": "shade_smooth", "blendername": "Shade Smooth", "type": bool, "default": True }
    ]
    style="ball_and_stick"
    # fmt: on

    def __init__(
        self,
        quality: int = 2,
        sphere_radii: float = 0.3,
        bond_split: bool = False,
        bond_find: bool = False,
        bond_radius: float = 0.3,
        color_blur: bool = False,
        shade_smooth: bool = True,
    ):
        self.quality = quality
        self.sphere_radii = sphere_radii
        self.bond_split = bond_split
        self.bond_find = bond_find
        self.bond_radius = bond_radius
        self.color_blur = color_blur
        self.shade_smooth = shade_smooth


class StyleCartoon(StyleBase):
    # fmt: off
    portdata: PortDataList = [
        { "name": "quality", "blendername": "Quality", "type": int, "default": 2 },
        { "name": "dssp", "blendername": "Peptide DSSP", "type": bool, "default": False },
        { "name": "cylinders", "blendername": "Peptide Cylinders", "type": bool, "default": False },
        { "name": "arrows", "blendername": "Peptide Arrows", "type": bool, "default": True },
        { "name": "rounded", "blendername": "Peptide Rounded", "type": bool, "default": False },
        { "name": "thickness", "blendername": "Peptide Thickness", "type": float, "default": 0.6000000238418579 },
        { "name": "width", "blendername": "Peptide Width", "type": float, "default": 2.200000047683716 },
        { "name": "loop_radius", "blendername": "Peptide Loop Radius", "type": float, "default": 0.30000001192092896 },
        { "name": "smoothing", "blendername": "Peptide Smoothing", "type": float, "default": 0.5 },
        { "name": "backbone_shape", "blendername": "Backbone Shape", "type": str, "default": "Cylinder" },
        { "name": "backbone_radius", "blendername": "Backbone Radius", "type": float, "default": 2.0 },
        { "name": "base_shape", "blendername": "Base Shape", "type": str, "default": "Rectangle" },
        { "name": "base_realize", "blendername": "Base Realize", "type": bool, "default": False },
        { "name": "color_blur", "blendername": "Color Blur", "type": bool, "default": True },
        { "name": "shade_smooth", "blendername": "Shade Smooth", "type": bool, "default": True }
    ]
    # fmt: on
    style="cartoon"

    def __init__(
        self,
        quality: int = 2,
        dssp: bool = False,
        cylinders: bool = False,
        arrows: bool = True,
        rounded: bool = False,
        thickness: float = 0.6000000238418579,
        width: float = 2.200000047683716,
        loop_radius: float = 0.30000001192092896,
        smoothing: float = 0.5,
        backbone_shape: str = "Cylinder",
        backbone_radius: float = 2.0,
        base_shape: str = "Rectangle",
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
    # fmt: off
    portdata: PortDataList = [
        {"name": "quality", "blendername": "Quality", "type": int, "default": 3},
        {"name": "radius", "blendername": "Radius", "type": float, "default": 1.6},
        {"name": "smoothing", "blendername": "Smoothing", "type": float, "default": 0.6},
        {"name": "color_blur", "blendername": "Color Blur", "type": bool, "default": False},
        {"name": "shade_smooth", "blendername": "Shade Smooth", "type": bool, "default": False},
        {"name": "backbone_smoothing", "blendername": "Backbone Smoothing", "type": float, "default": 0.5},
        {"name": "backbone_threshold", "blendername": "Backbone Threshold", "type": float, "default": 4.5},
        {"name": "backbone_radius", "blendername": "Backbone Radius", "type": float, "default": 1.6},
        {"name": "backbone_shape", "blendername": "Backbone Shape", "type": str, "default": "Cylinder"},
        {"name": "base_resolution", "blendername": "Base Resolution", "type": int, "default": 4},
        {"name": "base_realize", "blendername": "Base Realize", "type": bool, "default": False},
        {"name": "uv_map", "blendername": "UV Map", "type": bool, "default": False},
        {"name": "u_component_factor", "blendername": "U Component Factor", "type": float, "default": None},
    ]

    # fmt: on
    style = "ribbon"
    def __init__(
        self,
        quality: int = 3,
        radius: float = 1.6,
        smoothing: float = 0.6,
        color_blur: bool = False,
        shade_smooth: bool = False,
        backbone_smoothing: float = 0.5,
        backbone_threshold: float = 4.5,
        backbone_radius: float = 1.6,
        backbone_shape: str = "Cylinder",
        base_resolution: int = 4,
        base_realize: bool = False,
        uv_map: bool = False,
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
        self.backbone_shape = backbone_shape
        self.base_resolution = base_resolution
        self.base_realize = base_realize
        self.uv_map = uv_map
        self.u_component_factor = u_component_factor



class StyleSpheres(StyleBase):
    # fmt: off
    portdata: PortDataList = [
        {"name": "geometry", "blendername": "Sphere Geometry", "type": str, "default": "Point"},
        {"name": "radii", "blendername": "Sphere Radii", "type": float, "default": 0.800000011920929},
        {"name": "sphere_subdivisions", "blendername": "Sphere Subdivisions", "type": int, "default": 2},
        {"name": "shade_smooth", "blendername": "Shade Smooth", "type": bool, "default": True},
    ]
    # fmt: on
    style = "spheres"

    def __init__(
        self,
        geometry: str = "Point",  # make enum
        radii: float = 0.800000011920929,
        sphere_subdivisions: int = 2,
        shade_smooth: bool = True,
    ):
        self.geometry = geometry
        self.radii = radii
        self.sphere_subdivisions = sphere_subdivisions
        self.shade_smooth = shade_smooth


class StyleSticks(StyleBase):
    # fmt: off
    portdata: PortDataList = [
        { "name": "quality", "blendername": "Quality", "type": int, "default": 3 },
        { "name": "radius", "blendername": "Radius", "type": float, "default": 0.20000000298023224 },
        { "name": "color_blur", "blendername": "Color Blur", "type": bool, "default": False },
        { "name": "shade_smooth", "blendername": "Shade Smooth", "type": bool, "default": True }
    ]
    # fmt: on
    style = "sticks"

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
    # fmt: off
    portdata: PortDataList = [
        { "name": "quality", "blendername": "Quality", "type": int, "default": 3 },
        { "name": "scale_radii", "blendername": "Scale Radii", "type": float, "default": 1.5 },
        { "name": "probe_size", "blendername": "Probe Size", "type": float, "default": 1.0 },
        { "name": "relaxation_steps", "blendername": "Relaxation Steps", "type": int, "default": 10 },
        { "name": "separate", "blendername": "Separate By", "type": str, "default": "chain_id" },
        { "name": "group_id", "blendername": "Group ID", "type": int, "default": 0 },
        { "name": "color_source", "blendername": "Color Source", "type": str, "default": "Alpha Carbon" },
        { "name": "blur", "blendername": "Color Blur", "type": int, "default": 2 },
        { "name": "shade_smooth", "blendername": "Shade Smooth", "type": bool, "default": True }
    ]
    # fmt: on
    style = "surface"

    def __init__(
        self,
        quality: int = 3,
        scale_radii: float = 1.5,
        probe_size: float = 1.0,
        relaxation_steps: int = 10,
        separate: str = "chain_id",
        group_id: int = 0,
        color_source: str = "Alpha Carbon",
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
