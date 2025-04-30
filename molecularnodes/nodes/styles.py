from dataclasses import dataclass, replace, field, fields
from typing import List, Tuple, Union, Any, Dict

StyleClass = Union[
    "StyleBallandStick",
    "StyleCartoon",
    "StyleRibbon",
    "StyleSpheres",
    "StyleSticks",
    "StyleSurface",
]

# Define the type for a single port data entry for clarity
PortDataEntry = Dict[str, Any]
# Define the type for the list of port data entries
PortDataList = List[PortDataEntry]


class StyleBase:
    # Base class for styles. Subclasses will define their own portdata attribute.
    portdata: PortDataList = []


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
    portdata: PortDataList = [
        {"name": "quality", "blendername": "Quality", "type": int, "default": 3},
        {"name": "radius", "blendername": "Radius", "type": float, "default": 1.6},
        {
            "name": "smoothing",
            "blendername": "Smoothing",
            "type": float,
            "default": 0.6,
        },
        {
            "name": "color_blur",
            "blendername": "Color Blur",
            "type": bool,
            "default": False,
        },
        {
            "name": "shade_smooth",
            "blendername": "Shade Smooth",
            "type": bool,
            "default": False,
        },
    ]

    def __init__(
        self,
        style: str = "ribbon",
        quality: int = 3,
        radius: float = 1.6,
        smoothing: float = 0.6,
        color_blur: bool = False,
        shade_smooth: bool = False,
    ):
        self.style = style
        self.quality = quality
        self.radius = radius
        self.smoothing = smoothing
        self.color_blur = color_blur
        self.shade_smooth = shade_smooth


class StyleSpheres(StyleBase):
    portdata: PortDataList = [
        # fmt: off
        {
            "name": "geometry",
            "blendername": "Sphere Geometry",
            "type": bool,
            "default": "Point",
        },
        {"name": "radii", "blendername": "Sphere Radii", "type": float, "default": 0.8},
        {
            "name": "subdivisions",
            "blendername": "Sphere Subdivisions",
            "type": int,
            "default": 2,
        },
        {
            "name": "shade_smooth",
            "blendername": "Shade Smooth",
            "type": bool,
            "default": False,
        },
        # fmt: on
    ]

    def __init__(
        self,
        style: str = "sphere",
        as_mesh: bool = True,
        radii: float = 0.8,
        subdivisions: int = 2,
        shade_smooth: bool = False,
    ):
        self.style = style
        self.as_mesh = as_mesh  # This corresponds to a geometry option, but the portdata has a different key "geometry"
        self.radii = radii
        self.subdivisions = subdivisions
        self.shade_smooth = shade_smooth


class StyleSticks(StyleBase):
    # fmt: off
    portdata: PortDataList = [
        { "name": "quality", "blendername": "Quality", "type": int, "default": 2 },
        { "name": "radius", "blendername": "Radius", "type": float, "default": 0.2 },
        { "name": "color_blur", "blendername": "Color Blur", "type": bool, "default": False },
        { "name": "shade_smooth", "blendername": "Shade Smooth", "type": bool, "default": False }
    ]
    # fmt: on

    def __init__(
        self,
        style: str = "sticks",
        quality: int = 2,
        radius: float = 0.2,
        color_blur: bool = False,
        shade_smooth: bool = False,
    ):
        self.style = style
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

    def __init__(
        self,
        style: str = "surface",
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
        self.style = style
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
