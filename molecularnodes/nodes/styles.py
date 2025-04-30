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
        { "name": "as_mesh", "blendername": "As Mesh", "type": bool, "default": True },
        { "name": "sphere_radii", "blendername": "Sphere Radii", "type": float, "default": 0.3 },
        { "name": "bond_split", "blendername": "Bond Split", "type": bool, "default": False },
        { "name": "bond_find", "blendername": "Bond Find", "type": bool, "default": True },
        { "name": "bond_radius", "blendername": "Bond Radius", "type": float, "default": 0.3 },
        { "name": "color_blur", "blendername": "Color Blur", "type": bool, "default": False },
        { "name": "shade_smooth", "blendername": "Shade Smooth", "type": bool, "default": True }
    ]

    # Quality 2
    # Sphere Geometry Instance
    # Sphere Radii 0.30000001192092896
    # Bond Split False
    # Bond Find False
    # Bond Radius 0.30000001192092896
    # Color Blur False
    # Shade Smooth True
    style="ball_and_stick"
    # fmt: on

    def __init__(
        self,
        quality: int = 2,
        as_mesh: bool = True,
        sphere_radii: float = 0.3,
        bond_split: bool = False,
        bond_find: bool = True,
        bond_radius: float = 0.3,
        color_blur: bool = False,
        shade_smooth: bool = True,
    ):
        self.quality = quality
        self.as_mesh = as_mesh
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
        { "name": "dssp", "blendername": "DSSP", "type": bool, "default": False },
        { "name": "cylinders", "blendername": "Cylinders", "type": bool, "default": False },
        { "name": "arrows", "blendername": "Arrows", "type": bool, "default": True },
        { "name": "rounded", "blendername": "Rounded", "type": bool, "default": False },
        { "name": "thickness", "blendername": "Thickness", "type": float, "default": 0.6 },
        { "name": "width", "blendername": "Width", "type": float, "default": 2.2 },
        { "name": "loop_radius", "blendername": "Loop Radius", "type": float, "default": 0.3 },
        { "name": "smoothing", "blendername": "Smoothing", "type": float, "default": 0.5 },
        { "name": "color_blur", "blendername": "Color Blur", "type": bool, "default": False },
        { "name": "shade_smooth", "blendername": "Shade Smooth", "type": bool, "default": True }
    ]
    # fmt: on
    style="cartoon"

    # Quality 2
    # Peptide DSSP False
    # Peptide Cylinders False
    # Peptide Arrows True
    # Peptide Rounded False
    # Peptide Thickness 0.6000000238418579
    # Peptide Width 2.200000047683716
    # Peptide Loop Radius 0.30000001192092896
    # Peptide Smoothing 0.5
    # Backbone Shape Cylinder
    # Backbone Radius 2.0
    # Base Shape Rectangle
    # Base Realize False
    # Color Blur True
    # Shade Smooth True


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
        color_blur: bool = False,
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
    portdata: PortDataList = [
        # fmt: off
        {"name": "geometry", "blendername": "Sphere Geometry", "type": bool, "default": "Point"},
        {"name": "radii", "blendername": "Sphere Radii", "type": float, "default": 0.8},
        {"name": "sphere_subdivisions", "blendername": "Sphere Subdivisions", "type": int, "default": 2},
        {"name": "shade_smooth", "blendername": "Shade Smooth", "type": bool, "default": False},
        # fmt: on
    ]
    style = "spheres"


    # Sphere Geometry Point
    # Sphere Radii 0.800000011920929
    # Sphere Subdivisions 2
    # Shade Smooth True

    def __init__(
        self,
        geometry: str = "Point",  # make enum
        radii: float = 0.8,
        sphere_subdivisions: int = 2,
        shade_smooth: bool = False,
    ):
        self.geometry = geometry
        self.radii = radii
        self.sphere_subdivisions = sphere_subdivisions
        self.shade_smooth = shade_smooth


class StyleSticks(StyleBase):
    # fmt: off
    portdata: PortDataList = [
        { "name": "quality", "blendername": "Quality", "type": int, "default": 2 },
        { "name": "radius", "blendername": "Radius", "type": float, "default": 0.2 },
        { "name": "color_blur", "blendername": "Color Blur", "type": bool, "default": False },
        { "name": "shade_smooth", "blendername": "Shade Smooth", "type": bool, "default": False }
    ]
    style = "sticks"

    # Quality 3
    # Radius 0.20000000298023224
    # Color Blur False
    # Shade Smooth True

    # fmt: on
    def __init__(
        self,
        quality: int = 2,
        radius: float = 0.2,
        color_blur: bool = False,
        shade_smooth: bool = False,
    ):
        self.quality = quality
        self.radius = radius
        self.color_blur = color_blur
        self.shade_smooth = shade_smooth


class StyleSurface(StyleBase):
    # fmt: off
    portdata: PortDataList = [
        { "name": "quality", "blendername": "Quality", "type": int, "default": 3 },
        { "name": "separate", "blendername": "Separate", "type": bool, "default": True },
        { "name": "attribute", "blendername": "Attribute", "type": str, "default": "chain_id" },
        { "name": "scale_radii", "blendername": "Scale Radii", "type": float, "default": 1.5 },
        { "name": "probe_size", "blendername": "Probe Size", "type": float, "default": 1.0 },
        { "name": "triangulate", "blendername": "Triangulate", "type": bool, "default": False },
        { "name": "relaxation_steps", "blendername": "Relaxation Steps", "type": int, "default": 10 },
        { "name": "by_ca", "blendername": "by CA", "type": bool, "default": False },
        { "name": "blur", "blendername": "Blur", "type": int, "default": 2 },
        { "name": "shade_smooth", "blendername": "Shade Smooth", "type": bool, "default": True }
    ]
    # fmt: on
    style = "surface"

    # Quality 3
    # Scale Radii 1.5
    # Probe Size 1.0
    # Relaxation Steps 10
    # Separate By chain_id
    # Group ID 0
    # Color Source Alpha Carbon
    # Color Blur 2
    # Shade Smooth True
    def __init__(
        self,
        quality: int = 3,
        separate: bool = True,
        attribute: str = "chain_id",
        scale_radii: float = 1.5,
        probe_size: float = 1.0,
        triangulate: bool = False,
        relaxation_steps: int = 10,
        by_ca: bool = False,
        blur: int = 2,
        shade_smooth: bool = True,
    ):
        self.quality = quality
        self.separate = separate
        self.attribute = attribute
        self.scale_radii = scale_radii
        self.probe_size = probe_size
        self.triangulate = triangulate
        self.relaxation_steps = relaxation_steps
        self.by_ca = by_ca
        self.blur = blur
        self.shade_smooth = shade_smooth
