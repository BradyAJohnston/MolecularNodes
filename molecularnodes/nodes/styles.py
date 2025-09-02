"""
Auto-generated Style Classes

This file contains Style classes automatically generated from a Blender file.
Each class represents a Style node and provides a Python interface to configure the node parameters.

Generated classes follow the same pattern as the existing styles in molecularnodes.nodes.styles,
using the Socket dataclass system and StyleBase inheritance.
"""

from dataclasses import dataclass
from typing import Any, List, Optional, Tuple
from bpy.types import GeometryNodeGroup

__all__ = [
    "StyleBallAndStick",
    "StyleCartoon",
    "StyleDensitySurface",
    "StyleDensityWire",
    "StylePreset1",
    "StylePreset2",
    "StylePreset3",
    "StylePreset4",
    "StyleRibbon",
    "StyleSpheres",
    "StyleSticks",
    "StyleSurface",
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

        Parameters
        ----------
        node_style : GeometryNodeGroup
            The Blender GeometryNodeGroup to update
        """
        for input in node_style.inputs:
            if input.type != "GEOMETRY":
                for arg in self.socketdata:
                    if input.name == arg.blendername:
                        input.default_value = getattr(self, arg.name)


class StyleBallAndStick(StyleBase):
    """Style class for Style Ball and Stick

    Parameters
    ----------
    quality : int
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    sphere_geometry : Any
        Show spheres as a _Point Cloud_, _Instances_ of a mesh Icosphere, or realised _Mesh_ instances of an Icosphere. Point cloud is best for performance and should definitely be used if rendering in Cycles.
    sphere_radius : float
        Scale the `vdw_radii` attribute before setting the radius for the spheres
    bond_split : bool
        Split apart double and triple bonds visually
    bond_radius : float
        Set the radius for the generated bonds in Angstroms
    bond_find : bool
        Find possible bonds for the selected atoms based on a distance search. Unselected atoms maintain any bonds they already have. Bonds that are found are all treated as single bonds
    bond_find_scale : float
        Scale the VDW radii of the atoms when searching for bonds
    color_blur : bool
        Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
    shade_smooth : bool
        Apply smooth shading to the created geometry
    """

    style = "ball_and_stick"
    socketdata: SocketInfo = [
        Socket(name="quality", blendername="Quality"),
        Socket(name="sphere_geometry", blendername="Sphere Geometry"),
        Socket(name="sphere_radius", blendername="Sphere Radius"),
        Socket(name="bond_split", blendername="Bond Split"),
        Socket(name="bond_radius", blendername="Bond Radius"),
        Socket(name="bond_find", blendername="Bond Find"),
        Socket(name="bond_find_scale", blendername="Bond Find Scale"),
        Socket(name="color_blur", blendername="Color Blur"),
        Socket(name="shade_smooth", blendername="Shade Smooth"),
    ]

    def __init__(
        self,
        quality: int = 2,  # A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
        sphere_geometry: Any = "Instance",  # Show spheres as a _Point Cloud_, _Instances_ of a mesh Icosphere, or realised _Mesh_ instances of an Icosphere. Point cloud is best for performance and should definitely be used if rendering in Cycles.
        sphere_radius: float = 0.3,  # Scale the `vdw_radii` attribute before setting the radius for the spheres
        bond_split: bool = False,  # Split apart double and triple bonds visually
        bond_radius: float = 0.3,  # Set the radius for the generated bonds in Angstroms
        bond_find: bool = False,  # Find possible bonds for the selected atoms based on a distance search. Unselected atoms maintain any bonds they already have. Bonds that are found are all treated as single bonds
        bond_find_scale: float = 1.0,  # Scale the VDW radii of the atoms when searching for bonds
        color_blur: bool = False,  # Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
        shade_smooth: bool = True,  # Apply smooth shading to the created geometry
    ):
        """Style class for Style Ball and Stick

            Parameters
            ----------
        quality : int
            A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
        sphere_geometry : Any
            Show spheres as a _Point Cloud_, _Instances_ of a mesh Icosphere, or realised _Mesh_ instances of an Icosphere. Point cloud is best for performance and should definitely be used if rendering in Cycles.
        sphere_radius : float
            Scale the `vdw_radii` attribute before setting the radius for the spheres
        bond_split : bool
            Split apart double and triple bonds visually
        bond_radius : float
            Set the radius for the generated bonds in Angstroms
        bond_find : bool
            Find possible bonds for the selected atoms based on a distance search. Unselected atoms maintain any bonds they already have. Bonds that are found are all treated as single bonds
        bond_find_scale : float
            Scale the VDW radii of the atoms when searching for bonds
        color_blur : bool
            Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
        shade_smooth : bool
            Apply smooth shading to the created geometry
        """
        self.quality = quality
        self.sphere_geometry = sphere_geometry
        self.sphere_radius = sphere_radius
        self.bond_split = bond_split
        self.bond_radius = bond_radius
        self.bond_find = bond_find
        self.bond_find_scale = bond_find_scale
        self.color_blur = color_blur
        self.shade_smooth = shade_smooth


class StyleCartoon(StyleBase):
    """Style class for Style Cartoon

    Parameters
    ----------
    quality : int
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    peptide_dssp : bool
        Use the DSSP algorithm to compute the `sec_struct` attribute
    peptide_cylinders : bool
        Use cylinders for helices instead of ribbons
    peptide_arrows : bool
        User arrows for sheets
    peptide_rounded : bool
        Create rounded sheets and helices
    peptide_thickness : float
        Thickness for the sheets and helices
    peptide_width : float
        Width for the sheets and helices
    peptide_loop_radius : float
        Radius of the loops for unstructure regions
    peptide_smoothing : float
        Smoothing to apply to sheets
    backbone_shape : Any
        Value for Backbone Shape
    nucleic_width : float
        Value for Nucleic Width
    nucleic_thickness : float
        Value for Nucleic Thickness
    nucleic_radius : float
        Value for Nucleic Radius
    base_shape : Any
        Value for Base Shape
    base_realize : bool
        Value for Base Realize
    color_blur : bool
        Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
    shade_smooth : bool
        Apply smooth shading to the created geometry
    """

    style = "cartoon"
    socketdata: SocketInfo = [
        Socket(name="quality", blendername="Quality"),
        Socket(name="peptide_dssp", blendername="Peptide DSSP"),
        Socket(name="peptide_cylinders", blendername="Peptide Cylinders"),
        Socket(name="peptide_arrows", blendername="Peptide Arrows"),
        Socket(name="peptide_rounded", blendername="Peptide Rounded"),
        Socket(name="peptide_thickness", blendername="Peptide Thickness"),
        Socket(name="peptide_width", blendername="Peptide Width"),
        Socket(name="peptide_loop_radius", blendername="Peptide Loop Radius"),
        Socket(name="peptide_smoothing", blendername="Peptide Smoothing"),
        Socket(name="backbone_shape", blendername="Backbone Shape"),
        Socket(name="nucleic_width", blendername="Nucleic Width"),
        Socket(name="nucleic_thickness", blendername="Nucleic Thickness"),
        Socket(name="nucleic_radius", blendername="Nucleic Radius"),
        Socket(name="base_shape", blendername="Base Shape"),
        Socket(name="base_realize", blendername="Base Realize"),
        Socket(name="color_blur", blendername="Color Blur"),
        Socket(name="shade_smooth", blendername="Shade Smooth"),
    ]

    def __init__(
        self,
        quality: int = 2,  # A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
        peptide_dssp: bool = False,  # Use the DSSP algorithm to compute the `sec_struct` attribute
        peptide_cylinders: bool = False,  # Use cylinders for helices instead of ribbons
        peptide_arrows: bool = True,  # User arrows for sheets
        peptide_rounded: bool = False,  # Create rounded sheets and helices
        peptide_thickness: float = 0.6,  # Thickness for the sheets and helices
        peptide_width: float = 2.2,  # Width for the sheets and helices
        peptide_loop_radius: float = 0.3,  # Radius of the loops for unstructure regions
        peptide_smoothing: float = 0.5,  # Smoothing to apply to sheets
        backbone_shape: Any = "Cylinder",
        nucleic_width: float = 3.0,
        nucleic_thickness: float = 1.0,
        nucleic_radius: float = 2.0,
        base_shape: Any = "Rectangle",
        base_realize: bool = False,
        color_blur: bool = True,  # Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
        shade_smooth: bool = True,  # Apply smooth shading to the created geometry
    ):
        """Style class for Style Cartoon

            Parameters
            ----------
        quality : int
            A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
        peptide_dssp : bool
            Use the DSSP algorithm to compute the `sec_struct` attribute
        peptide_cylinders : bool
            Use cylinders for helices instead of ribbons
        peptide_arrows : bool
            User arrows for sheets
        peptide_rounded : bool
            Create rounded sheets and helices
        peptide_thickness : float
            Thickness for the sheets and helices
        peptide_width : float
            Width for the sheets and helices
        peptide_loop_radius : float
            Radius of the loops for unstructure regions
        peptide_smoothing : float
            Smoothing to apply to sheets
        backbone_shape : Any
            Value for Backbone Shape
        nucleic_width : float
            Value for Nucleic Width
        nucleic_thickness : float
            Value for Nucleic Thickness
        nucleic_radius : float
            Value for Nucleic Radius
        base_shape : Any
            Value for Base Shape
        base_realize : bool
            Value for Base Realize
        color_blur : bool
            Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
        shade_smooth : bool
            Apply smooth shading to the created geometry
        """
        self.quality = quality
        self.peptide_dssp = peptide_dssp
        self.peptide_cylinders = peptide_cylinders
        self.peptide_arrows = peptide_arrows
        self.peptide_rounded = peptide_rounded
        self.peptide_thickness = peptide_thickness
        self.peptide_width = peptide_width
        self.peptide_loop_radius = peptide_loop_radius
        self.peptide_smoothing = peptide_smoothing
        self.backbone_shape = backbone_shape
        self.nucleic_width = nucleic_width
        self.nucleic_thickness = nucleic_thickness
        self.nucleic_radius = nucleic_radius
        self.base_shape = base_shape
        self.base_realize = base_realize
        self.color_blur = color_blur
        self.shade_smooth = shade_smooth


class StyleDensitySurface(StyleBase):
    """Style class for Style Density Surface

    Parameters
    ----------
    threshold : float
        Value for Threshold
    shade_smooth : bool
        Apply smooth shading to the created geometry
    hide_dust : float
        Value for Hide Dust
    color : Tuple[float, float, float, float]
        Value for Color
    """

    style = "density_surface"
    socketdata: SocketInfo = [
        Socket(name="threshold", blendername="Threshold"),
        Socket(name="shade_smooth", blendername="Shade Smooth"),
        Socket(name="hide_dust", blendername="Hide Dust"),
        Socket(name="color", blendername="Color"),
    ]

    def __init__(
        self,
        threshold: float = 0.8,
        shade_smooth: bool = True,  # Apply smooth shading to the created geometry
        hide_dust: float = 0.0,
        color: Tuple[float, float, float, float] = (0.2, 0.51, 0.13, 1.0),
    ):
        """Style class for Style Density Surface

            Parameters
            ----------
        threshold : float
            Value for Threshold
        shade_smooth : bool
            Apply smooth shading to the created geometry
        hide_dust : float
            Value for Hide Dust
        color : Tuple[float, float, float, float]
            Value for Color
        """
        self.threshold = threshold
        self.shade_smooth = shade_smooth
        self.hide_dust = hide_dust
        self.color = color


class StyleDensityWire(StyleBase):
    """Style class for Style Density Wire

    Parameters
    ----------
    threshold : float
        Value for Threshold
    hide_dust : float
        Value for Hide Dust
    wire_radius : float
        Radius of the created wire (in relative nm)
    wire_resolution : int
        Value for Wire Resolution
    color : Tuple[float, float, float, float]
        Value for Color
    """

    style = "density_wire"
    socketdata: SocketInfo = [
        Socket(name="threshold", blendername="Threshold"),
        Socket(name="hide_dust", blendername="Hide Dust"),
        Socket(name="wire_radius", blendername="Wire Radius"),
        Socket(name="wire_resolution", blendername="Wire Resolution"),
        Socket(name="color", blendername="Color"),
    ]

    def __init__(
        self,
        threshold: float = 0.8,
        hide_dust: float = 20.0,
        wire_radius: float = 1.0,  # Radius of the created wire (in relative nm)
        wire_resolution: int = 3,
        color: Tuple[float, float, float, float] = (0.1, 0.39, 0.1, 1.0),
    ):
        """Style class for Style Density Wire

            Parameters
            ----------
        threshold : float
            Value for Threshold
        hide_dust : float
            Value for Hide Dust
        wire_radius : float
            Radius of the created wire (in relative nm)
        wire_resolution : int
            Value for Wire Resolution
        color : Tuple[float, float, float, float]
            Value for Color
        """
        self.threshold = threshold
        self.hide_dust = hide_dust
        self.wire_radius = wire_radius
        self.wire_resolution = wire_resolution
        self.color = color


class StylePreset1(StyleBase):
    """Style class for Style Preset 1

    Parameters
    ----------
    quality : int
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    color_blur : bool
        Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
    shade_smooth : bool
        Apply smooth shading to the created geometry
    """

    style = "preset_1"
    socketdata: SocketInfo = [
        Socket(name="quality", blendername="Quality"),
        Socket(name="color_blur", blendername="Color Blur"),
        Socket(name="shade_smooth", blendername="Shade Smooth"),
    ]

    def __init__(
        self,
        quality: int = 3,  # A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
        color_blur: bool = False,  # Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
        shade_smooth: bool = True,  # Apply smooth shading to the created geometry
    ):
        """Style class for Style Preset 1

            Parameters
            ----------
        quality : int
            A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
        color_blur : bool
            Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
        shade_smooth : bool
            Apply smooth shading to the created geometry
        """
        self.quality = quality
        self.color_blur = color_blur
        self.shade_smooth = shade_smooth


class StylePreset2(StyleBase):
    """Style class for Style Preset 2

    Parameters
    ----------
    quality : int
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    color_blur : bool
        Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
    shade_smooth : bool
        Apply smooth shading to the created geometry
    """

    style = "preset_2"
    socketdata: SocketInfo = [
        Socket(name="quality", blendername="Quality"),
        Socket(name="color_blur", blendername="Color Blur"),
        Socket(name="shade_smooth", blendername="Shade Smooth"),
    ]

    def __init__(
        self,
        quality: int = 3,  # A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
        color_blur: bool = False,  # Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
        shade_smooth: bool = True,  # Apply smooth shading to the created geometry
    ):
        """Style class for Style Preset 2

            Parameters
            ----------
        quality : int
            A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
        color_blur : bool
            Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
        shade_smooth : bool
            Apply smooth shading to the created geometry
        """
        self.quality = quality
        self.color_blur = color_blur
        self.shade_smooth = shade_smooth


class StylePreset3(StyleBase):
    """Style class for Style Preset 3

    Parameters
    ----------
    quality : int
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    shade_smooth : bool
        Apply smooth shading to the created geometry
    """

    style = "preset_3"
    socketdata: SocketInfo = [
        Socket(name="quality", blendername="Quality"),
        Socket(name="shade_smooth", blendername="Shade Smooth"),
    ]

    def __init__(
        self,
        quality: int = 3,  # A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
        shade_smooth: bool = True,  # Apply smooth shading to the created geometry
    ):
        """Style class for Style Preset 3

            Parameters
            ----------
        quality : int
            A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
        shade_smooth : bool
            Apply smooth shading to the created geometry
        """
        self.quality = quality
        self.shade_smooth = shade_smooth


class StylePreset4(StyleBase):
    """Style class for Style Preset 4

    Parameters
    ----------
    quality : int
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    color_blur : bool
        Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
    shade_smooth : bool
        Apply smooth shading to the created geometry
    """

    style = "preset_4"
    socketdata: SocketInfo = [
        Socket(name="quality", blendername="Quality"),
        Socket(name="color_blur", blendername="Color Blur"),
        Socket(name="shade_smooth", blendername="Shade Smooth"),
    ]

    def __init__(
        self,
        quality: int = 3,  # A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
        color_blur: bool = True,  # Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
        shade_smooth: bool = True,  # Apply smooth shading to the created geometry
    ):
        """Style class for Style Preset 4

            Parameters
            ----------
        quality : int
            A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
        color_blur : bool
            Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
        shade_smooth : bool
            Apply smooth shading to the created geometry
        """
        self.quality = quality
        self.color_blur = color_blur
        self.shade_smooth = shade_smooth


class StyleRibbon(StyleBase):
    """Style class for Style Ribbon

    Parameters
    ----------
    quality : int
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    backbone_smoothing : float
        Smoothen the sheet ribbons such as beta-sheets
    backbone_threshold : float
        Distance (Angstroms) over which subsequent CA points are treated as a new chain
    backbone_radius : float
        Value for Backbone Radius
    nucleic_backbone_shape : Any
        Value for Nucleic Backbone Shape
    nucleic_backbone_radius : float
        Value for Nucleic Backbone Radius
    backbone_width : float
        Value for Backbone Width
    backbone_thickness : float
        Value for Backbone Thickness
    base_scale : Tuple[float, float, float]
        Value for Base Scale
    base_resolution : int
        Value for Base Resolution
    base_realize : bool
        Value for Base Realize
    uv_map : bool
        Compute and store the `uv_map` for the final protein ribbon geometry
    u_component : Any
        Value for U Component
    color_blur : bool
        Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
    shade_smooth : bool
        Apply smooth shading to the created geometry
    """

    style = "ribbon"
    socketdata: SocketInfo = [
        Socket(name="quality", blendername="Quality"),
        Socket(name="backbone_smoothing", blendername="Backbone Smoothing"),
        Socket(name="backbone_threshold", blendername="Backbone Threshold"),
        Socket(name="backbone_radius", blendername="Backbone Radius"),
        Socket(name="nucleic_backbone_shape", blendername="Nucleic Backbone Shape"),
        Socket(name="nucleic_backbone_radius", blendername="Nucleic Backbone Radius"),
        Socket(name="backbone_width", blendername="Backbone Width"),
        Socket(name="backbone_thickness", blendername="Backbone Thickness"),
        Socket(name="base_scale", blendername="Base Scale"),
        Socket(name="base_resolution", blendername="Base Resolution"),
        Socket(name="base_realize", blendername="Base Realize"),
        Socket(name="uv_map", blendername="UV Map"),
        Socket(name="u_component", blendername="U Component"),
        Socket(name="color_blur", blendername="Color Blur"),
        Socket(name="shade_smooth", blendername="Shade Smooth"),
    ]

    def __init__(
        self,
        quality: int = 3,  # A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
        backbone_smoothing: float = 0.5,  # Smoothen the sheet ribbons such as beta-sheets
        backbone_threshold: float = 4.5,  # Distance (Angstroms) over which subsequent CA points are treated as a new chain
        backbone_radius: float = 1.6,
        nucleic_backbone_shape: Any = "Cylinder",
        nucleic_backbone_radius: float = 1.6,
        backbone_width: float = 3.0,
        backbone_thickness: float = 1.0,
        base_scale: Tuple[float, float, float] = (2.5, 0.5, 7.0),
        base_resolution: int = 4,
        base_realize: bool = False,
        uv_map: bool = False,  # Compute and store the `uv_map` for the final protein ribbon geometry
        u_component: Any = "Factor",
        color_blur: bool = False,  # Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
        shade_smooth: bool = True,  # Apply smooth shading to the created geometry
    ):
        """Style class for Style Ribbon

            Parameters
            ----------
        quality : int
            A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
        backbone_smoothing : float
            Smoothen the sheet ribbons such as beta-sheets
        backbone_threshold : float
            Distance (Angstroms) over which subsequent CA points are treated as a new chain
        backbone_radius : float
            Value for Backbone Radius
        nucleic_backbone_shape : Any
            Value for Nucleic Backbone Shape
        nucleic_backbone_radius : float
            Value for Nucleic Backbone Radius
        backbone_width : float
            Value for Backbone Width
        backbone_thickness : float
            Value for Backbone Thickness
        base_scale : Tuple[float, float, float]
            Value for Base Scale
        base_resolution : int
            Value for Base Resolution
        base_realize : bool
            Value for Base Realize
        uv_map : bool
            Compute and store the `uv_map` for the final protein ribbon geometry
        u_component : Any
            Value for U Component
        color_blur : bool
            Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
        shade_smooth : bool
            Apply smooth shading to the created geometry
        """
        self.quality = quality
        self.backbone_smoothing = backbone_smoothing
        self.backbone_threshold = backbone_threshold
        self.backbone_radius = backbone_radius
        self.nucleic_backbone_shape = nucleic_backbone_shape
        self.nucleic_backbone_radius = nucleic_backbone_radius
        self.backbone_width = backbone_width
        self.backbone_thickness = backbone_thickness
        self.base_scale = base_scale
        self.base_resolution = base_resolution
        self.base_realize = base_realize
        self.uv_map = uv_map
        self.u_component = u_component
        self.color_blur = color_blur
        self.shade_smooth = shade_smooth


class StyleSpheres(StyleBase):
    """Style class for Style Spheres

    Parameters
    ----------
    geometry : Any
        Show spheres as a _Point Cloud_, _Instances_ of a mesh Icosphere, or realised _Mesh_ instances of an Icosphere. Point cloud is best for performance and should definitely be used if rendering in Cycles.
    radius : float
        Scale the `vdw_radii` of the atom when setting the radius of the spheres
    subdivisions : int
        Number of subdicisions when using _Instances_ or _Mesh_ to represent atoms
    shade_smooth : bool
        Apply smooth shading when using _Instances_ or _Mesh_
    """

    style = "spheres"
    socketdata: SocketInfo = [
        Socket(name="geometry", blendername="Geometry"),
        Socket(name="radius", blendername="Radius"),
        Socket(name="subdivisions", blendername="Subdivisions"),
        Socket(name="shade_smooth", blendername="Shade Smooth"),
    ]

    def __init__(
        self,
        geometry: Any = "Point",  # Show spheres as a _Point Cloud_, _Instances_ of a mesh Icosphere, or realised _Mesh_ instances of an Icosphere. Point cloud is best for performance and should definitely be used if rendering in Cycles.
        radius: float = 0.8,  # Scale the `vdw_radii` of the atom when setting the radius of the spheres
        subdivisions: int = 2,  # Number of subdicisions when using _Instances_ or _Mesh_ to represent atoms
        shade_smooth: bool = True,  # Apply smooth shading when using _Instances_ or _Mesh_
    ):
        """Style class for Style Spheres

            Parameters
            ----------
        geometry : Any
            Show spheres as a _Point Cloud_, _Instances_ of a mesh Icosphere, or realised _Mesh_ instances of an Icosphere. Point cloud is best for performance and should definitely be used if rendering in Cycles.
        radius : float
            Scale the `vdw_radii` of the atom when setting the radius of the spheres
        subdivisions : int
            Number of subdicisions when using _Instances_ or _Mesh_ to represent atoms
        shade_smooth : bool
            Apply smooth shading when using _Instances_ or _Mesh_
        """
        self.geometry = geometry
        self.radius = radius
        self.subdivisions = subdivisions
        self.shade_smooth = shade_smooth


class StyleSticks(StyleBase):
    """Style class for Style Sticks

    Parameters
    ----------
    quality : int
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    radius : float
        Radius of the sticks in Angstroms
    color_blur : bool
        Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
    shade_smooth : bool
        Apply smooth shading to the created geometry
    """

    style = "sticks"
    socketdata: SocketInfo = [
        Socket(name="quality", blendername="Quality"),
        Socket(name="radius", blendername="Radius"),
        Socket(name="color_blur", blendername="Color Blur"),
        Socket(name="shade_smooth", blendername="Shade Smooth"),
    ]

    def __init__(
        self,
        quality: int = 3,  # A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
        radius: float = 0.2,  # Radius of the sticks in Angstroms
        color_blur: bool = False,  # Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
        shade_smooth: bool = True,  # Apply smooth shading to the created geometry
    ):
        """Style class for Style Sticks

            Parameters
            ----------
        quality : int
            A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
        radius : float
            Radius of the sticks in Angstroms
        color_blur : bool
            Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
        shade_smooth : bool
            Apply smooth shading to the created geometry
        """
        self.quality = quality
        self.radius = radius
        self.color_blur = color_blur
        self.shade_smooth = shade_smooth


class StyleSurface(StyleBase):
    """Style class for Style Surface

    Parameters
    ----------
    quality : int
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    scale_radius : float
        Scale the VDW radii of the atoms when creating the surface
    probe_size : float
        Size of the probe that is used to check for solvent accessibility (Angstroms)
    relaxation_steps : int
        Number of times smoothening is applied to the generate surface stretched between the atoms
    separate_by : Any
        Value for Separate By
    group_id : int
        Value for Group ID
    color_source : Any
        Value for Color Source
    color_blur : int
        Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
    shade_smooth : bool
        Apply smooth shading to the created geometry
    """

    style = "surface"
    socketdata: SocketInfo = [
        Socket(name="quality", blendername="Quality"),
        Socket(name="scale_radius", blendername="Scale Radius"),
        Socket(name="probe_size", blendername="Probe Size"),
        Socket(name="relaxation_steps", blendername="Relaxation Steps"),
        Socket(name="separate_by", blendername="Separate By"),
        Socket(name="group_id", blendername="Group ID"),
        Socket(name="color_source", blendername="Color Source"),
        Socket(name="color_blur", blendername="Color Blur"),
        Socket(name="shade_smooth", blendername="Shade Smooth"),
    ]

    def __init__(
        self,
        quality: int = 3,  # A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
        scale_radius: float = 1.5,  # Scale the VDW radii of the atoms when creating the surface
        probe_size: float = 1.0,  # Size of the probe that is used to check for solvent accessibility (Angstroms)
        relaxation_steps: int = 10,  # Number of times smoothening is applied to the generate surface stretched between the atoms
        separate_by: Any = "chain_id",
        group_id: int = 0,
        color_source: Any = "Alpha Carbon",
        color_blur: int = 2,  # Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
        shade_smooth: bool = True,  # Apply smooth shading to the created geometry
    ):
        """Style class for Style Surface

            Parameters
            ----------
        quality : int
            A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
        scale_radius : float
            Scale the VDW radii of the atoms when creating the surface
        probe_size : float
            Size of the probe that is used to check for solvent accessibility (Angstroms)
        relaxation_steps : int
            Number of times smoothening is applied to the generate surface stretched between the atoms
        separate_by : Any
            Value for Separate By
        group_id : int
            Value for Group ID
        color_source : Any
            Value for Color Source
        color_blur : int
            Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
        shade_smooth : bool
            Apply smooth shading to the created geometry
        """
        self.quality = quality
        self.scale_radius = scale_radius
        self.probe_size = probe_size
        self.relaxation_steps = relaxation_steps
        self.separate_by = separate_by
        self.group_id = group_id
        self.color_source = color_source
        self.color_blur = color_blur
        self.shade_smooth = shade_smooth
