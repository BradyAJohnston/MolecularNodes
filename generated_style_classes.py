"""
Auto-generated Style Classes

This file contains Style classes automatically generated from the MN_data_file_4.4.blend file.
Each class represents a Style node and provides a Python interface to configure the node parameters.

Generated classes follow the same pattern as the existing styles in molecularnodes.nodes.styles,
using the Socket dataclass system and StyleBase inheritance.
"""

from dataclasses import dataclass
from enum import Enum
from typing import List, Optional, Tuple, Any
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

        Args:
            node_style: The Blender GeometryNodeGroup to update
        """
        for input in node_style.inputs:
            if input.type != "GEOMETRY":
                for arg in self.socketdata:
                    if input.name == arg.blendername:
                        input.default_value = getattr(self, arg.name)


class SphereGeometryEnum(str, Enum):
    """Show spheres as a _Point Cloud_, _Instances_ of a mesh Icosphere, or realised _Mesh_ instances of an Icosphere. Point cloud is best for performance and should definitely be used if rendering in Cycles.
    
    Options
    -------
    Point : Point - Point cloud representation
    Instance : Instance - Instanced mesh spheres
    Mesh : Mesh - Realized mesh spheres
    """
    POINT = "Point"
    INSTANCE = "Instance"
    MESH = "Mesh"
class StyleBallAndStick(StyleBase):
    """Style class for Style Ball and Stick
    
    Parameters
    ----------
    quality : int
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    selection : bool
        Selection of atoms to apply this style to, discarding unselected points
    sphere_geometry : SphereGeometryEnum
        Show spheres as a _Point Cloud_, _Instances_ of a mesh Icosphere, or realised _Mesh_ instances of an Icosphere. Point cloud is best for performance and should definitely be used if rendering in Cycles.. Options: 'Point', 'Instance', 'Mesh'
    sphere_radii : float
        Scale the `vdw_radii` attribute before setting the radius for the spheres
    bond_split : bool
        Split apart double and triple bonds visually
    bond_find : bool
        Find possible bonds for the selected atoms based on a distance search. Unselected atoms maintain any bonds they already have. Bonds that are found are all treated as single bonds
    bond_radius : float
        Set the radius for the generated bonds in Angstroms
    color_blur : bool
        Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
    shade_smooth : bool
        Apply smooth shading to the created geometry
    material : Any
        Material to apply to the resulting geometry
    """
    style = "ball_and_stick"
    socketdata: SocketInfo = [
        Socket(name="quality", blendername="Quality"),
        Socket(name="selection", blendername="Selection"),
        Socket(name="sphere_geometry", blendername="Sphere Geometry"),
        Socket(name="sphere_radii", blendername="Sphere Radii"),
        Socket(name="bond_split", blendername="Bond Split"),
        Socket(name="bond_find", blendername="Bond Find"),
        Socket(name="bond_radius", blendername="Bond Radius"),
        Socket(name="color_blur", blendername="Color Blur"),
        Socket(name="shade_smooth", blendername="Shade Smooth"),
        Socket(name="material", blendername="Material"),
    ]

    def __init__(
        self,
        quality: int = 2,  # A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
        selection: bool = True,  # Selection of atoms to apply this style to, discarding unselected points
        sphere_geometry: SphereGeometryEnum = 'Instance',  # Show spheres as a _Point Cloud_, _Instances_ of a mesh Icosphere, or realised _Mesh_ instances of an Icosphere. Point cloud is best for performance and should definitely be used if rendering in Cycles.
        sphere_radii: float = 0.30000001192092896,  # Scale the `vdw_radii` attribute before setting the radius for the spheres
        bond_split: bool = False,  # Split apart double and triple bonds visually
        bond_find: bool = False,  # Find possible bonds for the selected atoms based on a distance search. Unselected atoms maintain any bonds they already have. Bonds that are found are all treated as single bonds
        bond_radius: float = 0.30000001192092896,  # Set the radius for the generated bonds in Angstroms
        color_blur: bool = False,  # Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
        shade_smooth: bool = True,  # Apply smooth shading to the created geometry
        material: Any = None,  # Material to apply to the resulting geometry
    ):
        """Style class for Style Ball and Stick

        Parameters
        ----------
    quality : int
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    selection : bool
        Selection of atoms to apply this style to, discarding unselected points
    sphere_geometry : SphereGeometryEnum
        Show spheres as a _Point Cloud_, _Instances_ of a mesh Icosphere, or realised _Mesh_ instances of an Icosphere. Point cloud is best for performance and should definitely be used if rendering in Cycles.. Options: 'Point', 'Instance', 'Mesh'
    sphere_radii : float
        Scale the `vdw_radii` attribute before setting the radius for the spheres
    bond_split : bool
        Split apart double and triple bonds visually
    bond_find : bool
        Find possible bonds for the selected atoms based on a distance search. Unselected atoms maintain any bonds they already have. Bonds that are found are all treated as single bonds
    bond_radius : float
        Set the radius for the generated bonds in Angstroms
    color_blur : bool
        Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
    shade_smooth : bool
        Apply smooth shading to the created geometry
    material : Any
        Material to apply to the resulting geometry
        """
        self.quality = quality
        self.selection = selection
        self.sphere_geometry = sphere_geometry
        self.sphere_radii = sphere_radii
        self.bond_split = bond_split
        self.bond_find = bond_find
        self.bond_radius = bond_radius
        self.color_blur = color_blur
        self.shade_smooth = shade_smooth
        self.material = material


class BackboneShapeEnum(str, Enum):
    """Enum for Backbone Shape in Style Cartoon
    
    Options
    -------
    Cylinder : Cylinder - Cylindrical backbone
    Rectangle : Rectangle - Rectangular backbone
    """
    CYLINDER = "Cylinder"
    RECTANGLE = "Rectangle"
class BaseShapeEnum(str, Enum):
    """Enum for Base Shape in Style Cartoon
    
    Options
    -------
    Rectangle : Rectangle - Rectangular base shape
    Cylinder : Cylinder - Cylindrical base shape
    Square : Square - Square base shape
    """
    RECTANGLE = "Rectangle"
    CYLINDER = "Cylinder"
    SQUARE = "Square"
class StyleCartoon(StyleBase):
    """Style class for Style Cartoon
    
    Parameters
    ----------
    selection : bool
        Selection of atoms to apply this style to, discarding unselected points
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
    backbone_shape : BackboneShapeEnum
        Value for Backbone Shape. Options: 'Cylinder', 'Rectangle'
    backbone_radius : float
        Value for Backbone Radius
    base_shape : BaseShapeEnum
        Value for Base Shape. Options: 'Rectangle', 'Cylinder', 'Square'
    base_realize : bool
        Value for Base Realize
    color_blur : bool
        Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
    shade_smooth : bool
        Apply smooth shading to the created geometry
    material : Any
        Material to apply to the resulting geometry
    """
    style = "cartoon"
    socketdata: SocketInfo = [
        Socket(name="selection", blendername="Selection"),
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
        Socket(name="backbone_radius", blendername="Backbone Radius"),
        Socket(name="base_shape", blendername="Base Shape"),
        Socket(name="base_realize", blendername="Base Realize"),
        Socket(name="color_blur", blendername="Color Blur"),
        Socket(name="shade_smooth", blendername="Shade Smooth"),
        Socket(name="material", blendername="Material"),
    ]

    def __init__(
        self,
        selection: bool = True,  # Selection of atoms to apply this style to, discarding unselected points
        quality: int = 2,  # A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
        peptide_dssp: bool = False,  # Use the DSSP algorithm to compute the `sec_struct` attribute
        peptide_cylinders: bool = False,  # Use cylinders for helices instead of ribbons
        peptide_arrows: bool = True,  # User arrows for sheets
        peptide_rounded: bool = False,  # Create rounded sheets and helices
        peptide_thickness: float = 0.6000000238418579,  # Thickness for the sheets and helices
        peptide_width: float = 2.200000047683716,  # Width for the sheets and helices
        peptide_loop_radius: float = 0.30000001192092896,  # Radius of the loops for unstructure regions
        peptide_smoothing: float = 0.5,  # Smoothing to apply to sheets
        backbone_shape: BackboneShapeEnum = 'Cylinder',
        backbone_radius: float = 2.0,
        base_shape: BaseShapeEnum = 'Rectangle',
        base_realize: bool = False,
        color_blur: bool = True,  # Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
        shade_smooth: bool = True,  # Apply smooth shading to the created geometry
        material: Any = None,  # Material to apply to the resulting geometry
    ):
        """Style class for Style Cartoon

        Parameters
        ----------
    selection : bool
        Selection of atoms to apply this style to, discarding unselected points
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
    backbone_shape : BackboneShapeEnum
        Value for Backbone Shape. Options: 'Cylinder', 'Rectangle'
    backbone_radius : float
        Value for Backbone Radius
    base_shape : BaseShapeEnum
        Value for Base Shape. Options: 'Rectangle', 'Cylinder', 'Square'
    base_realize : bool
        Value for Base Realize
    color_blur : bool
        Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
    shade_smooth : bool
        Apply smooth shading to the created geometry
    material : Any
        Material to apply to the resulting geometry
        """
        self.selection = selection
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
        self.backbone_radius = backbone_radius
        self.base_shape = base_shape
        self.base_realize = base_realize
        self.color_blur = color_blur
        self.shade_smooth = shade_smooth
        self.material = material


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
    material : Any
        Material to apply to the resulting geometry
    """
    style = "density_surface"
    socketdata: SocketInfo = [
        Socket(name="threshold", blendername="Threshold"),
        Socket(name="shade_smooth", blendername="Shade Smooth"),
        Socket(name="hide_dust", blendername="Hide Dust"),
        Socket(name="color", blendername="Color"),
        Socket(name="material", blendername="Material"),
    ]

    def __init__(
        self,
        threshold: float = 0.800000011920929,
        shade_smooth: bool = True,  # Apply smooth shading to the created geometry
        hide_dust: float = 0.0,
        color: Tuple[float, float, float, float] = (0.1994359940290451, 0.5091630220413208, 0.13218000531196594, 1.0),
        material: Any = None,  # Material to apply to the resulting geometry
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
    material : Any
        Material to apply to the resulting geometry
        """
        self.threshold = threshold
        self.shade_smooth = shade_smooth
        self.hide_dust = hide_dust
        self.color = color
        self.material = material


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
    material : Any
        Material to apply to the resulting geometry
    """
    style = "density_wire"
    socketdata: SocketInfo = [
        Socket(name="threshold", blendername="Threshold"),
        Socket(name="hide_dust", blendername="Hide Dust"),
        Socket(name="wire_radius", blendername="Wire Radius"),
        Socket(name="wire_resolution", blendername="Wire Resolution"),
        Socket(name="color", blendername="Color"),
        Socket(name="material", blendername="Material"),
    ]

    def __init__(
        self,
        threshold: float = 0.800000011920929,
        hide_dust: float = 20.0,
        wire_radius: float = 1.0,  # Radius of the created wire (in relative nm)
        wire_resolution: int = 3,
        color: Tuple[float, float, float, float] = (0.10174982994794846, 0.3931145668029785, 0.10474135726690292, 1.0),
        material: Any = None,  # Material to apply to the resulting geometry
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
    material : Any
        Material to apply to the resulting geometry
        """
        self.threshold = threshold
        self.hide_dust = hide_dust
        self.wire_radius = wire_radius
        self.wire_resolution = wire_resolution
        self.color = color
        self.material = material


class StylePreset1(StyleBase):
    """Style class for Style Preset 1
    
    Parameters
    ----------
    selection : bool
        Selection of atoms to apply this style to, discarding unselected points
    quality : int
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    color_blur : bool
        Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
    shade_smooth : bool
        Apply smooth shading to the created geometry
    material : Any
        Material to apply to the resulting geometry
    """
    style = "preset_1"
    socketdata: SocketInfo = [
        Socket(name="selection", blendername="Selection"),
        Socket(name="quality", blendername="Quality"),
        Socket(name="color_blur", blendername="Color Blur"),
        Socket(name="shade_smooth", blendername="Shade Smooth"),
        Socket(name="material", blendername="Material"),
    ]

    def __init__(
        self,
        selection: bool = True,  # Selection of atoms to apply this style to, discarding unselected points
        quality: int = 3,  # A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
        color_blur: bool = False,  # Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
        shade_smooth: bool = True,  # Apply smooth shading to the created geometry
        material: Any = None,  # Material to apply to the resulting geometry
    ):
        """Style class for Style Preset 1

        Parameters
        ----------
    selection : bool
        Selection of atoms to apply this style to, discarding unselected points
    quality : int
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    color_blur : bool
        Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
    shade_smooth : bool
        Apply smooth shading to the created geometry
    material : Any
        Material to apply to the resulting geometry
        """
        self.selection = selection
        self.quality = quality
        self.color_blur = color_blur
        self.shade_smooth = shade_smooth
        self.material = material


class StylePreset2(StyleBase):
    """Style class for Style Preset 2
    
    Parameters
    ----------
    selection : bool
        Selection of atoms to apply this style to, discarding unselected points
    quality : int
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    color_blur : bool
        Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
    shade_smooth : bool
        Apply smooth shading to the created geometry
    material : Any
        Material to apply to the resulting geometry
    """
    style = "preset_2"
    socketdata: SocketInfo = [
        Socket(name="selection", blendername="Selection"),
        Socket(name="quality", blendername="Quality"),
        Socket(name="color_blur", blendername="Color Blur"),
        Socket(name="shade_smooth", blendername="Shade Smooth"),
        Socket(name="material", blendername="Material"),
    ]

    def __init__(
        self,
        selection: bool = True,  # Selection of atoms to apply this style to, discarding unselected points
        quality: int = 3,  # A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
        color_blur: bool = False,  # Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
        shade_smooth: bool = True,  # Apply smooth shading to the created geometry
        material: Any = None,  # Material to apply to the resulting geometry
    ):
        """Style class for Style Preset 2

        Parameters
        ----------
    selection : bool
        Selection of atoms to apply this style to, discarding unselected points
    quality : int
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    color_blur : bool
        Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
    shade_smooth : bool
        Apply smooth shading to the created geometry
    material : Any
        Material to apply to the resulting geometry
        """
        self.selection = selection
        self.quality = quality
        self.color_blur = color_blur
        self.shade_smooth = shade_smooth
        self.material = material


class StylePreset3(StyleBase):
    """Style class for Style Preset 3
    
    Parameters
    ----------
    selection : bool
        Selection of atoms to apply this style to, discarding unselected points
    quality : int
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    shade_smooth : bool
        Apply smooth shading to the created geometry
    material : Any
        Material to apply to the resulting geometry
    """
    style = "preset_3"
    socketdata: SocketInfo = [
        Socket(name="selection", blendername="Selection"),
        Socket(name="quality", blendername="Quality"),
        Socket(name="shade_smooth", blendername="Shade Smooth"),
        Socket(name="material", blendername="Material"),
    ]

    def __init__(
        self,
        selection: bool = True,  # Selection of atoms to apply this style to, discarding unselected points
        quality: int = 3,  # A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
        shade_smooth: bool = True,  # Apply smooth shading to the created geometry
        material: Any = None,  # Material to apply to the resulting geometry
    ):
        """Style class for Style Preset 3

        Parameters
        ----------
    selection : bool
        Selection of atoms to apply this style to, discarding unselected points
    quality : int
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    shade_smooth : bool
        Apply smooth shading to the created geometry
    material : Any
        Material to apply to the resulting geometry
        """
        self.selection = selection
        self.quality = quality
        self.shade_smooth = shade_smooth
        self.material = material


class StylePreset4(StyleBase):
    """Style class for Style Preset 4
    
    Parameters
    ----------
    selection : bool
        Selection of atoms to apply this style to, discarding unselected points
    quality : int
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    color_blur : bool
        Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
    shade_smooth : bool
        Apply smooth shading to the created geometry
    material : Any
        Material to apply to the resulting geometry
    """
    style = "preset_4"
    socketdata: SocketInfo = [
        Socket(name="selection", blendername="Selection"),
        Socket(name="quality", blendername="Quality"),
        Socket(name="color_blur", blendername="Color Blur"),
        Socket(name="shade_smooth", blendername="Shade Smooth"),
        Socket(name="material", blendername="Material"),
    ]

    def __init__(
        self,
        selection: bool = True,  # Selection of atoms to apply this style to, discarding unselected points
        quality: int = 3,  # A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
        color_blur: bool = True,  # Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
        shade_smooth: bool = True,  # Apply smooth shading to the created geometry
        material: Any = None,  # Material to apply to the resulting geometry
    ):
        """Style class for Style Preset 4

        Parameters
        ----------
    selection : bool
        Selection of atoms to apply this style to, discarding unselected points
    quality : int
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    color_blur : bool
        Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
    shade_smooth : bool
        Apply smooth shading to the created geometry
    material : Any
        Material to apply to the resulting geometry
        """
        self.selection = selection
        self.quality = quality
        self.color_blur = color_blur
        self.shade_smooth = shade_smooth
        self.material = material


class BackboneShapeEnum(str, Enum):
    """Enum for Backbone Shape in Style Ribbon
    
    Options
    -------
    Cylinder : Cylinder - Cylindrical backbone
    Rectangle : Rectangle - Rectangular backbone
    """
    CYLINDER = "Cylinder"
    RECTANGLE = "Rectangle"
class UComponentEnum(str, Enum):
    """Enum for U Component in Style Ribbon
    
    Options
    -------
    Factor : Factor - Use factor component
    Position : Position - Use position component
    """
    FACTOR = "Factor"
    POSITION = "Position"
class StyleRibbon(StyleBase):
    """Style class for Style Ribbon
    
    Parameters
    ----------
    selection : bool
        Selection of atoms to apply this style to, discarding unselected points
    quality : int
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    backbone_smoothing : float
        Smoothen the sheet ribbons such as beta-sheets
    backbone_threshold : float
        Distance (Angstroms) over which subsequent CA points are treated as a new chain
    backbone_radius : float
        Value for Backbone Radius
    backbone_shape : BackboneShapeEnum
        Value for Backbone Shape. Options: 'Cylinder', 'Rectangle'
    backbone_radius : float
        Value for Backbone Radius
    base_scale : Tuple[float, float, float]
        Value for Base Scale
    base_resolution : int
        Value for Base Resolution
    base_realize : bool
        Value for Base Realize
    uv_map : bool
        Compute and store the `uv_map` for the final protein ribbon geometry
    u_component : UComponentEnum
        Value for U Component. Options: 'Factor', 'Position'
    color_blur : bool
        Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
    shade_smooth : bool
        Apply smooth shading to the created geometry
    material : Any
        Material to apply to the resulting geometry
    """
    style = "ribbon"
    socketdata: SocketInfo = [
        Socket(name="selection", blendername="Selection"),
        Socket(name="quality", blendername="Quality"),
        Socket(name="backbone_smoothing", blendername="Backbone Smoothing"),
        Socket(name="backbone_threshold", blendername="Backbone Threshold"),
        Socket(name="backbone_radius", blendername="Backbone Radius"),
        Socket(name="backbone_shape", blendername="Backbone Shape"),
        Socket(name="backbone_radius", blendername="Backbone Radius"),
        Socket(name="base_scale", blendername="Base Scale"),
        Socket(name="base_resolution", blendername="Base Resolution"),
        Socket(name="base_realize", blendername="Base Realize"),
        Socket(name="uv_map", blendername="UV Map"),
        Socket(name="u_component", blendername="U Component"),
        Socket(name="color_blur", blendername="Color Blur"),
        Socket(name="shade_smooth", blendername="Shade Smooth"),
        Socket(name="material", blendername="Material"),
    ]

    def __init__(
        self,
        selection: bool = True,  # Selection of atoms to apply this style to, discarding unselected points
        quality: int = 3,  # A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
        backbone_smoothing: float = 0.5,  # Smoothen the sheet ribbons such as beta-sheets
        backbone_threshold: float = 4.5,  # Distance (Angstroms) over which subsequent CA points are treated as a new chain
        backbone_radius: float = 1.600000023841858,
        backbone_shape: BackboneShapeEnum = 'Cylinder',
        backbone_radius: float = 2.0,
        base_scale: Tuple[float, float, float] = (2.5, 0.5, 7.0),
        base_resolution: int = 4,
        base_realize: bool = False,
        uv_map: bool = False,  # Compute and store the `uv_map` for the final protein ribbon geometry
        u_component: UComponentEnum = 'Factor',
        color_blur: bool = True,  # Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
        shade_smooth: bool = True,  # Apply smooth shading to the created geometry
        material: Any = None,  # Material to apply to the resulting geometry
    ):
        """Style class for Style Ribbon

        Parameters
        ----------
    selection : bool
        Selection of atoms to apply this style to, discarding unselected points
    quality : int
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    backbone_smoothing : float
        Smoothen the sheet ribbons such as beta-sheets
    backbone_threshold : float
        Distance (Angstroms) over which subsequent CA points are treated as a new chain
    backbone_radius : float
        Value for Backbone Radius
    backbone_shape : BackboneShapeEnum
        Value for Backbone Shape. Options: 'Cylinder', 'Rectangle'
    backbone_radius : float
        Value for Backbone Radius
    base_scale : Tuple[float, float, float]
        Value for Base Scale
    base_resolution : int
        Value for Base Resolution
    base_realize : bool
        Value for Base Realize
    uv_map : bool
        Compute and store the `uv_map` for the final protein ribbon geometry
    u_component : UComponentEnum
        Value for U Component. Options: 'Factor', 'Position'
    color_blur : bool
        Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
    shade_smooth : bool
        Apply smooth shading to the created geometry
    material : Any
        Material to apply to the resulting geometry
        """
        self.selection = selection
        self.quality = quality
        self.backbone_smoothing = backbone_smoothing
        self.backbone_threshold = backbone_threshold
        self.backbone_radius = backbone_radius
        self.backbone_shape = backbone_shape
        self.backbone_radius = backbone_radius
        self.base_scale = base_scale
        self.base_resolution = base_resolution
        self.base_realize = base_realize
        self.uv_map = uv_map
        self.u_component = u_component
        self.color_blur = color_blur
        self.shade_smooth = shade_smooth
        self.material = material


class SphereGeometryEnum(str, Enum):
    """Show spheres as a _Point Cloud_, _Instances_ of a mesh Icosphere, or realised _Mesh_ instances of an Icosphere. Point cloud is best for performance and should definitely be used if rendering in Cycles.
    
    Options
    -------
    Point : Point - Point cloud representation
    Instance : Instance - Instanced mesh spheres
    Mesh : Mesh - Realized mesh spheres
    """
    POINT = "Point"
    INSTANCE = "Instance"
    MESH = "Mesh"
class StyleSpheres(StyleBase):
    """Style class for Style Spheres
    
    Parameters
    ----------
    selection : bool
        Selection of atoms to apply this style to, discarding unselected points
    sphere_geometry : SphereGeometryEnum
        Show spheres as a _Point Cloud_, _Instances_ of a mesh Icosphere, or realised _Mesh_ instances of an Icosphere. Point cloud is best for performance and should definitely be used if rendering in Cycles.. Options: 'Point', 'Instance', 'Mesh'
    sphere_radii : float
        Scale the `vdw_radii` of the atom when setting the radius of the spheres
    sphere_subdivisions : int
        Number of subdicisions when using _Instances_ or _Mesh_ to represent atoms
    shade_smooth : bool
        Apply smooth shading when using _Instances_ or _Mesh_
    material : Any
        Material to apply to the resulting geometry
    """
    style = "spheres"
    socketdata: SocketInfo = [
        Socket(name="selection", blendername="Selection"),
        Socket(name="sphere_geometry", blendername="Sphere Geometry"),
        Socket(name="sphere_radii", blendername="Sphere Radii"),
        Socket(name="sphere_subdivisions", blendername="Sphere Subdivisions"),
        Socket(name="shade_smooth", blendername="Shade Smooth"),
        Socket(name="material", blendername="Material"),
    ]

    def __init__(
        self,
        selection: bool = True,  # Selection of atoms to apply this style to, discarding unselected points
        sphere_geometry: SphereGeometryEnum = 'Point',  # Show spheres as a _Point Cloud_, _Instances_ of a mesh Icosphere, or realised _Mesh_ instances of an Icosphere. Point cloud is best for performance and should definitely be used if rendering in Cycles.
        sphere_radii: float = 0.800000011920929,  # Scale the `vdw_radii` of the atom when setting the radius of the spheres
        sphere_subdivisions: int = 2,  # Number of subdicisions when using _Instances_ or _Mesh_ to represent atoms
        shade_smooth: bool = True,  # Apply smooth shading when using _Instances_ or _Mesh_
        material: Any = None,  # Material to apply to the resulting geometry
    ):
        """Style class for Style Spheres

        Parameters
        ----------
    selection : bool
        Selection of atoms to apply this style to, discarding unselected points
    sphere_geometry : SphereGeometryEnum
        Show spheres as a _Point Cloud_, _Instances_ of a mesh Icosphere, or realised _Mesh_ instances of an Icosphere. Point cloud is best for performance and should definitely be used if rendering in Cycles.. Options: 'Point', 'Instance', 'Mesh'
    sphere_radii : float
        Scale the `vdw_radii` of the atom when setting the radius of the spheres
    sphere_subdivisions : int
        Number of subdicisions when using _Instances_ or _Mesh_ to represent atoms
    shade_smooth : bool
        Apply smooth shading when using _Instances_ or _Mesh_
    material : Any
        Material to apply to the resulting geometry
        """
        self.selection = selection
        self.sphere_geometry = sphere_geometry
        self.sphere_radii = sphere_radii
        self.sphere_subdivisions = sphere_subdivisions
        self.shade_smooth = shade_smooth
        self.material = material


class StyleSticks(StyleBase):
    """Style class for Style Sticks
    
    Parameters
    ----------
    selection : bool
        Selection of atoms to apply this style to, discarding unselected points
    quality : int
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    radius : float
        Radius of the sticks in Angstroms
    color_blur : bool
        Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
    shade_smooth : bool
        Apply smooth shading to the created geometry
    material : Any
        Material to apply to the resulting geometry
    """
    style = "sticks"
    socketdata: SocketInfo = [
        Socket(name="selection", blendername="Selection"),
        Socket(name="quality", blendername="Quality"),
        Socket(name="radius", blendername="Radius"),
        Socket(name="color_blur", blendername="Color Blur"),
        Socket(name="shade_smooth", blendername="Shade Smooth"),
        Socket(name="material", blendername="Material"),
    ]

    def __init__(
        self,
        selection: bool = True,  # Selection of atoms to apply this style to, discarding unselected points
        quality: int = 3,  # A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
        radius: float = 0.20000000298023224,  # Radius of the sticks in Angstroms
        color_blur: bool = False,  # Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
        shade_smooth: bool = True,  # Apply smooth shading to the created geometry
        material: Any = None,  # Material to apply to the resulting geometry
    ):
        """Style class for Style Sticks

        Parameters
        ----------
    selection : bool
        Selection of atoms to apply this style to, discarding unselected points
    quality : int
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    radius : float
        Radius of the sticks in Angstroms
    color_blur : bool
        Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
    shade_smooth : bool
        Apply smooth shading to the created geometry
    material : Any
        Material to apply to the resulting geometry
        """
        self.selection = selection
        self.quality = quality
        self.radius = radius
        self.color_blur = color_blur
        self.shade_smooth = shade_smooth
        self.material = material


class SeparateByEnum(str, Enum):
    """Enum for Separate By in Style Surface
    
    Options
    -------
    chain_id : Chain ID - Separate by chain identifier
    res_id : Residue ID - Separate by residue identifier
    entity_id : Entity ID - Separate by entity identifier
    """
    CHAIN_ID = "chain_id"
    RES_ID = "res_id"
    ENTITY_ID = "entity_id"
class ColorSourceEnum(str, Enum):
    """Enum for Color Source in Style Surface
    
    Options
    -------
    Alpha Carbon : Alpha Carbon - Color from alpha carbon atoms
    Backbone : Backbone - Color from backbone atoms
    All : All - Color from all atoms
    """
    ALPHA_CARBON = "Alpha Carbon"
    BACKBONE = "Backbone"
    ALL = "All"
class StyleSurface(StyleBase):
    """Style class for Style Surface
    
    Parameters
    ----------
    selection : bool
        Selection of atoms to apply this style to, discarding unselected points
    quality : int
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    scale_radii : float
        Scale the VDW radii of the atoms when creating the surface
    probe_size : float
        Size of the probe that is used to check for solvent accessibility (Angstroms)
    relaxation_steps : int
        Number of times smoothening is applied to the generate surface stretched between the atoms
    separate_by : SeparateByEnum
        Value for Separate By. Options: 'chain_id', 'res_id', 'entity_id'
    group_id : int
        Value for Group ID
    color_source : ColorSourceEnum
        Value for Color Source. Options: 'Alpha Carbon', 'Backbone', 'All'
    color_blur : int
        Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
    shade_smooth : bool
        Apply smooth shading to the created geometry
    material : Any
        Material to apply to the resulting geometry
    """
    style = "surface"
    socketdata: SocketInfo = [
        Socket(name="selection", blendername="Selection"),
        Socket(name="quality", blendername="Quality"),
        Socket(name="scale_radii", blendername="Scale Radii"),
        Socket(name="probe_size", blendername="Probe Size"),
        Socket(name="relaxation_steps", blendername="Relaxation Steps"),
        Socket(name="separate_by", blendername="Separate By"),
        Socket(name="group_id", blendername="Group ID"),
        Socket(name="color_source", blendername="Color Source"),
        Socket(name="color_blur", blendername="Color Blur"),
        Socket(name="shade_smooth", blendername="Shade Smooth"),
        Socket(name="material", blendername="Material"),
    ]

    def __init__(
        self,
        selection: bool = True,  # Selection of atoms to apply this style to, discarding unselected points
        quality: int = 3,  # A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
        scale_radii: float = 1.5,  # Scale the VDW radii of the atoms when creating the surface
        probe_size: float = 1.0,  # Size of the probe that is used to check for solvent accessibility (Angstroms)
        relaxation_steps: int = 10,  # Number of times smoothening is applied to the generate surface stretched between the atoms
        separate_by: SeparateByEnum = 'chain_id',
        group_id: int = 0,
        color_source: ColorSourceEnum = 'Alpha Carbon',
        color_blur: int = 2,  # Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
        shade_smooth: bool = True,  # Apply smooth shading to the created geometry
        material: Any = None,  # Material to apply to the resulting geometry
    ):
        """Style class for Style Surface

        Parameters
        ----------
    selection : bool
        Selection of atoms to apply this style to, discarding unselected points
    quality : int
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    scale_radii : float
        Scale the VDW radii of the atoms when creating the surface
    probe_size : float
        Size of the probe that is used to check for solvent accessibility (Angstroms)
    relaxation_steps : int
        Number of times smoothening is applied to the generate surface stretched between the atoms
    separate_by : SeparateByEnum
        Value for Separate By. Options: 'chain_id', 'res_id', 'entity_id'
    group_id : int
        Value for Group ID
    color_source : ColorSourceEnum
        Value for Color Source. Options: 'Alpha Carbon', 'Backbone', 'All'
    color_blur : int
        Interpolate between colors when enabled. When disabled the faces will take their color from their corresponding atom without interpolating
    shade_smooth : bool
        Apply smooth shading to the created geometry
    material : Any
        Material to apply to the resulting geometry
        """
        self.selection = selection
        self.quality = quality
        self.scale_radii = scale_radii
        self.probe_size = probe_size
        self.relaxation_steps = relaxation_steps
        self.separate_by = separate_by
        self.group_id = group_id
        self.color_source = color_source
        self.color_blur = color_blur
        self.shade_smooth = shade_smooth
        self.material = material
