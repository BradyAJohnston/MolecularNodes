"""
Style Classes

MolecularNodes uses Geometry Nodes to define molecule representations within blender. However,
it is desireable to make these nodes accessible for scripting. The classes defined in this file
provide a way to script against the GeometryNode nodetrees. When invoked on a style these classes
will simply override the base styles.

Because the nodetrees in the MN_data_file_{version}.blend are the source of truth, we add
function to parity between the nodetrees and the class representations.

"""

from dataclasses import dataclass
from typing import Any, ClassVar, Dict, List, Optional, Tuple
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

class SocketInfo(List[Dict[str, str]]):
    pass

class StyleBase:
    style: ClassVar[str]
    portdata: ClassVar[List[Dict[str, str]]]

    def update_style_node(self, node_style: GeometryNodeGroup):
        for input in node_style.inputs:
            if input.type != "GEOMETRY":
                for arg in self.portdata:
                    name = arg["name"]
                    blendername = arg.get("blendername")
                    if input.name == blendername:
                        input.default_value = getattr(self, name)

    def get_blender_name(self, prop_name: str) -> str:
        """Get the Blender name for a given property name."""
        for port in self.portdata:
            if port["name"] == prop_name:
                return port["blendername"]
        return prop_name  # Return the original name if not found

    @classmethod
    def get_name_mapping(cls) -> Dict[str, str]:
        """Returns a dictionary mapping property names to Blender names."""
        return {port["name"]: port["blendername"] for port in cls.portdata}


@dataclass
class StyleBallandStick(StyleBase):
    quality: int = 2
    geometry: str = "Instance"  # enum
    sphere_radii: float = 0.3
    bond_split: bool = False
    bond_find: bool = False
    bond_radius: float = 0.3
    color_blur: bool = False
    shade_smooth: bool = True

    # Class variables
    style: ClassVar[str] = "ball_and_stick"
    portdata: ClassVar[List[Dict[str, str]]] = [
        {"name": "quality", "blendername": "Quality"},
        {"name": "geometry", "blendername": "Sphere Geometry"},
        {"name": "sphere_radii", "blendername": "Sphere Radii"},
        {"name": "bond_split", "blendername": "Bond Split"},
        {"name": "bond_find", "blendername": "Bond Find"},
        {"name": "bond_radius", "blendername": "Bond Radius"},
        {"name": "color_blur", "blendername": "Color Blur"},
        {"name": "shade_smooth", "blendername": "Shade Smooth"},
    ]


@dataclass
class StyleCartoon(StyleBase):
    quality: int = 2
    dssp: bool = False
    cylinders: bool = False
    arrows: bool = True
    rounded: bool = False
    thickness: float = 0.6
    width: float = 2.2
    loop_radius: float = 0.3
    smoothing: float = 0.5
    backbone_shape: str = "Cylinder"  # enum
    backbone_radius: float = 2.0
    base_shape: str = "Rectangle"  # enum
    base_realize: bool = False
    color_blur: bool = True
    shade_smooth: bool = True

    # Class variables
    style: ClassVar[str] = "cartoon"
    portdata: ClassVar[List[Dict[str, str]]] = [
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

@dataclass
class StyleRibbon(StyleBase):
    quality: int = 3
    color_blur: bool = False
    shade_smooth: bool = False
    backbone_smoothing: float = 0.5
    backbone_threshold: float = 4.5
    base_scale: Tuple[float, float, float] = (2.5, 0.5, 7.0)
    backbone_radius: float = 2.0
    backbone_shape: str = "Cylinder"
    base_resolution: int = 4
    base_realize: bool = False
    uv_map: bool = False
    u_component: Optional[Any] = None  # Using Any as generic type

    # Class variables
    style: ClassVar[str] = "ribbon"
    portdata: ClassVar[List[Dict[str, str]]] = [
        {"name": "quality", "blendername": "Quality"},
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
    ]

@dataclass
class StyleSpheres(StyleBase):
    geometry: str = "Point"  # enum: "Point" (Point Cloud), "Instances" (Instances of a mesh Icosphere), or "Mesh" (realised Mesh)
    radii: float = 0.8
    sphere_subdivisions: int = 2
    shade_smooth: bool = True

    # Class variables
    style: ClassVar[str] = "spheres"
    portdata: ClassVar[List[Dict[str, str]]] = [
        {"name": "geometry", "blendername": "Sphere Geometry"},
        {"name": "radii", "blendername": "Sphere Radii"},
        {"name": "sphere_subdivisions", "blendername": "Sphere Subdivisions"},
        {"name": "shade_smooth", "blendername": "Shade Smooth"},
    ]

@dataclass
class StyleSticks(StyleBase):
    quality: int = 3
    radius: float = 0.2
    color_blur: bool = False
    shade_smooth: bool = True

    # Class variables
    style: ClassVar[str] = "sticks"
    portdata: ClassVar[List[Dict[str, str]]] = [
        {"name": "quality", "blendername": "Quality"},
        {"name": "radius", "blendername": "Radius"},
        {"name": "color_blur", "blendername": "Color Blur"},
        {"name": "shade_smooth", "blendername": "Shade Smooth"},
    ]

@dataclass
class StyleSurface(StyleBase):
    quality: int = 3
    scale_radii: float = 1.5
    probe_size: float = 1.0
    relaxation_steps: int = 10
    separate: str = "chain_id"  # enum
    group_id: int = 0
    color_source: str = "Alpha Carbon"  # enum
    blur: int = 2
    shade_smooth: bool = True

    # Class variables
    style: ClassVar[str] = "surface"
    portdata: ClassVar[List[Dict[str, str]]] = [
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
