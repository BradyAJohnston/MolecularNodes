"""
Style Classes

MolecularNodes uses Geometry Nodes to define molecule representations within blender. However,
it is desireable to make these nodes accessible for scripting. The classes defined in this file
provide a way to script against the GeometryNode nodetrees. When invoked on a style these classes
will simply override the base styles.

Because the nodetrees in the MN_data_file_{version}.blend are the source of truth, we add
function to parity between the nodetrees and the class representations.

Each style uses a Socket dataclass to define the mapping between class attributes and
Blender node inputs, making it easier to maintain and extend the styles. The Socket
dataclass establishes a clear relationship between Python attribute names and the corresponding
Blender node input socket names, providing type safety and better IDE support than the
previous dictionary-based approach.
"""

from dataclasses import dataclass
from typing import List, Optional, Tuple
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


@dataclass
class Socket:
    """Represents a mapping between a class attribute and a Blender node input socket.

    Attributes:
        name: The name of the attribute in the Style class
        blendername: The corresponding name of the input socket in the Blender node
    """

    name: str
    blendername: Optional[str] = None


SocketInfo = List[Socket]


class StyleBase:
    style: str
    socketdata: SocketInfo = []

    def update_style_node(self, node_style: GeometryNodeGroup):
        """Update the Blender node inputs with values from this style's attributes.

        Args:
            node_style: The Blender GeometryNodeGroup to update
        """
        for input in node_style.inputs:
            if input.type != "GEOMETRY":
                for arg in self.socketdata:
                    if input.name == arg.blendername:
                        input.default_value = getattr(self, arg.name)


class StyleBallandStick(StyleBase):
    style = "ball_and_stick"
    socketdata: SocketInfo = [
        Socket(name="quality", blendername="Quality"),
        Socket(name="geometry", blendername="Sphere Geometry"),
        Socket(name="sphere_radii", blendername="Sphere Radii"),
        Socket(name="bond_split", blendername="Bond Split"),
        Socket(name="bond_find", blendername="Bond Find"),
        Socket(name="bond_radius", blendername="Bond Radius"),
        Socket(name="color_blur", blendername="Color Blur"),
        Socket(name="shade_smooth", blendername="Shade Smooth"),
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
    socketdata: SocketInfo = [
        Socket(name="quality", blendername="Quality"),
        Socket(name="dssp", blendername="Peptide DSSP"),
        Socket(name="cylinders", blendername="Peptide Cylinders"),
        Socket(name="arrows", blendername="Peptide Arrows"),
        Socket(name="rounded", blendername="Peptide Rounded"),
        Socket(name="thickness", blendername="Peptide Thickness"),
        Socket(name="width", blendername="Peptide Width"),
        Socket(name="loop_radius", blendername="Peptide Loop Radius"),
        Socket(name="smoothing", blendername="Peptide Smoothing"),
        Socket(name="backbone_shape", blendername="Backbone Shape"),
        Socket(name="backbone_radius", blendername="Backbone Radius"),
        Socket(name="base_shape", blendername="Base Shape"),
        Socket(name="base_realize", blendername="Base Realize"),
        Socket(name="color_blur", blendername="Color Blur"),
        Socket(name="shade_smooth", blendername="Shade Smooth"),
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
    socketdata: SocketInfo = [
        Socket(name="quality", blendername="Quality"),
        Socket(name="color_blur", blendername="Color Blur"),
        Socket(name="shade_smooth", blendername="Shade Smooth"),
        Socket(name="backbone_smoothing", blendername="Backbone Smoothing"),
        Socket(name="backbone_threshold", blendername="Backbone Threshold"),
        Socket(name="backbone_radius", blendername="Backbone Radius"),
        Socket(name="backbone_shape", blendername="Backbone Shape"),
        Socket(name="base_scale", blendername="Base Scale"),
        Socket(name="base_resolution", blendername="Base Resolution"),
        Socket(name="base_realize", blendername="Base Realize"),
        Socket(name="uv_map", blendername="UV Map"),
        Socket(name="u_component", blendername="U Component"),
    ]

    def __init__(
        self,
        quality: int = 3,
        color_blur: bool = False,
        shade_smooth: bool = False,
        backbone_smoothing: float = 0.5,
        backbone_threshold: float = 4.5,
        base_scale: Tuple[float, float, float] = (2.5, 0.5, 7.0),
        backbone_radius: float = 2.0,
        backbone_shape: str = "Cylinder",
        base_resolution: int = 4,
        base_realize: bool = False,
        uv_map: bool = False,
        u_component=None,
        # u_component_factor=None,
    ):
        self.quality = quality
        self.color_blur = color_blur
        self.shade_smooth = shade_smooth
        self.backbone_smoothing = backbone_smoothing
        self.backbone_threshold = backbone_threshold
        self.backbone_radius = backbone_radius
        self.backbone_shape = backbone_shape  # enum?
        self.base_scale = base_scale
        self.base_resolution = base_resolution
        self.base_realize = base_realize
        self.uv_map = uv_map
        self.u_component = u_component


class StyleSpheres(StyleBase):
    style = "spheres"
    socketdata: SocketInfo = [
        Socket(name="geometry", blendername="Sphere Geometry"),
        Socket(name="radii", blendername="Sphere Radii"),
        Socket(name="sphere_subdivisions", blendername="Sphere Subdivisions"),
        Socket(name="shade_smooth", blendername="Shade Smooth"),
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
    socketdata: SocketInfo = [
        Socket(name="quality", blendername="Quality"),
        Socket(name="radius", blendername="Radius"),
        Socket(name="color_blur", blendername="Color Blur"),
        Socket(name="shade_smooth", blendername="Shade Smooth"),
    ]

    def __init__(
        self,
        quality: int = 3,
        radius: float = 0.2,
        color_blur: bool = False,
        shade_smooth: bool = True,
    ):
        self.quality = quality
        self.radius = radius
        self.color_blur = color_blur
        self.shade_smooth = shade_smooth


class StyleSurface(StyleBase):
    style = "surface"
    socketdata: SocketInfo = [
        Socket(name="quality", blendername="Quality"),
        Socket(name="scale_radii", blendername="Scale Radii"),
        Socket(name="probe_size", blendername="Probe Size"),
        Socket(name="relaxation_steps", blendername="Relaxation Steps"),
        Socket(name="separate", blendername="Separate By"),
        Socket(name="group_id", blendername="Group ID"),
        Socket(name="color_source", blendername="Color Source"),
        Socket(name="blur", blendername="Color Blur"),
        Socket(name="shade_smooth", blendername="Shade Smooth"),
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
