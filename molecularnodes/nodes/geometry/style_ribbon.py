# Node-group asset 'Style Ribbon' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    CustomGeometryGroup,
    FloatSocket,
    GeometrySocket,
    IntegerSocket,
    MaterialSocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import (
    InputBoolean,
    InputFloat,
    InputGeometry,
    InputInteger,
    InputMaterial,
    InputMenu,
    InputVector,
)
from ._shared.cleanup import Cleanup
from ._shared.mn_utils_style_ribbon_nucleic import MN_utils_style_ribbon_nucleic
from .atoms_to_ca_curves import AtomsToCACurves
from .check_geometry import CheckGeometry
from .curve_custom_profile import CurveCustomProfile
from .curve_rotation import CurveRotation
from .evaluate_on_atoms import EvaluateOnAtoms
from .n2_index_angle import Group2IndexAngle
from .offset_index import OffsetIndex
from .separate_polymers import SeparatePolymers


class MN_utils_style_ribbon_peptide(CustomGeometryGroup):
    _name = ".MN_utils_style_ribbon_peptide"
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry._mn_utils_style_ribbon_peptide"}

    def _build_group(self, tree):
        atoms = tree.inputs.geometry(
            "Atoms", description="Atomic geometry that contains vertices and edges"
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Selection of atoms to apply this node to",
            hide_value=True,
        )
        quality = tree.inputs.integer("Quality", 3, min_value=0, max_value=6)
        bs_smoothing = tree.inputs.float(
            "BS Smoothing", 0.5, min_value=0.0, max_value=1.0, subtype="FACTOR"
        )
        radius = tree.inputs.float("Radius", 1.6, min_value=0.0)
        shade_smooth = tree.inputs.boolean(
            "Shade Smooth",
            True,
            description="Apply smooth shading to the created geometry",
        )
        material = tree.inputs.material(
            "Material", description="Material to apply to the resulting geometry"
        )
        uv_map = tree.inputs.boolean("UV Map", False)
        u_component = tree.inputs.menu("U Component", optional_label=True)
        threshold = tree.inputs.float(
            "Threshold", 4.5, min_value=0.0, max_value=10_000.0
        )
        geometry = tree.outputs.geometry("Geometry")
        curve = tree.outputs.geometry("Curve")

        group = AtomsToCACurves(
            atoms=atoms,
            selection=selection,
            bs_smoothing=bs_smoothing,
            threshold=threshold,
        )
        map_range = Group2IndexAngle(
            index_a=OffsetIndex(offset=-1), index_c=OffsetIndex(offset=1)
        ).o.angle.map_range(1.6899998, 2.861593, to_max=0.25)
        set_curve_radius = g.SetCurveRadius(
            curve=g.SetCurveNormal(curve=group), radius=radius
        )
        fillet_curve = g.FilletCurve(curve=set_curve_radius, radius=map_range)
        fillet_curve.node.mute = True
        group_1 = CurveCustomProfile(
            curve=g.SetHandleType(curve=g.SetSplineType.bezier(fillet_curve)),
            subdivisions=quality * g.Integer(integer=2),
            profile_type="Default Profile",
            uv_map=uv_map,
            u_component=u_component,
            socket_6=CurveRotation(),
            profile_resolution=quality * g.Integer(integer=3),
            input_14=0.0,
        )
        (
            Cleanup(
                geometry=group_1,
                color_source=group,
                material=material,
                shade_smooth=shade_smooth,
            )
            >> geometry
        )

        set_curve_radius >> curve

        u_component.default_value = "Factor"


class StyleRibbon(AssetGeometryGroup):
    """
    Style Ribbon

    Parameters
    ----------
    atoms : InputGeometry
        Atomic geometry that contains vertices and edges
    selection : InputBoolean
        Selection of atoms to apply this style to, discarding unselected points
    quality : InputInteger
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    peptide_radius : InputFloat
        Peptide Radius
    backbone_smoothing : InputFloat
        Smoothen the sheet ribbons such as beta-sheets
    backbone_threshold : InputFloat
        Distance (Angstroms) over which subsequent CA points are treated as a new chain
    uv_map : InputBoolean
        Compute and store the `uv_map` for the final protein ribbon geometry
    u_component : InputMenu | Literal["Factor", "Length"]
        Store either the 'Length' or the 'Factor' of the curve as the U component.
    nucleic_backbone_shape : InputMenu | Literal["Cicular", "Rectangular"]
        Nucleic Backbone Shape
    nucleic_backbone_radius : InputFloat
        Nucleic Backbone Radius
    nucleic_backbone_width : InputFloat
        Nucleic Backbone Width
    nucleic_backbone_thickness : InputFloat
        Nucleic Backbone Thickness
    base_geometry : InputGeometry
        Base Geometry
    base_scale : InputVector
        Base Scale
    base_resolution : InputInteger
        Base Resolution
    base_realize : InputBoolean
        Base Realize
    shade_smooth : InputBoolean
        Apply smooth shading to the created geometry
    material : InputMaterial
        Material to apply to the resulting geometry

    Inputs
    ------
    i.atoms : GeometrySocket
        Atomic geometry that contains vertices and edges
    i.selection : BooleanSocket
        Selection of atoms to apply this style to, discarding unselected points
    i.quality : IntegerSocket
        A lower value results in less geometry, with a higher value meaning better looking but more dense geometry
    i.peptide_radius : FloatSocket
        Peptide Radius
    i.backbone_smoothing : FloatSocket
        Smoothen the sheet ribbons such as beta-sheets
    i.backbone_threshold : FloatSocket
        Distance (Angstroms) over which subsequent CA points are treated as a new chain
    i.uv_map : BooleanSocket
        Compute and store the `uv_map` for the final protein ribbon geometry
    i.u_component : MenuSocket
        Store either the 'Length' or the 'Factor' of the curve as the U component.
    i.nucleic_backbone_shape : MenuSocket
        Nucleic Backbone Shape
    i.nucleic_backbone_radius : FloatSocket
        Nucleic Backbone Radius
    i.nucleic_backbone_width : FloatSocket
        Nucleic Backbone Width
    i.nucleic_backbone_thickness : FloatSocket
        Nucleic Backbone Thickness
    i.base_geometry : GeometrySocket
        Base Geometry
    i.base_scale : VectorSocket
        Base Scale
    i.base_resolution : IntegerSocket
        Base Resolution
    i.base_realize : BooleanSocket
        Base Realize
    i.shade_smooth : BooleanSocket
        Apply smooth shading to the created geometry
    i.material : MaterialSocket
        Material to apply to the resulting geometry

    Outputs
    -------
    o.geometry : GeometrySocket
        The generated geometry for the style node group
    """

    _name = "Style Ribbon"
    _asset_name = "Style Ribbon"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "node_tool_idname": "geometry.style_ribbon",
        "show_modifier_manage_panel": False,
        "is_modifier": True,
    }

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic geometry that contains vertices and edges"""
        selection: BooleanSocket
        """Selection of atoms to apply this style to, discarding unselected points"""
        quality: IntegerSocket
        """A lower value results in less geometry, with a higher value meaning better looking but more dense geometry"""
        peptide_radius: FloatSocket
        """Peptide Radius"""
        backbone_smoothing: FloatSocket
        """Smoothen the sheet ribbons such as beta-sheets"""
        backbone_threshold: FloatSocket
        """Distance (Angstroms) over which subsequent CA points are treated as a new chain"""
        uv_map: BooleanSocket
        """Compute and store the `uv_map` for the final protein ribbon geometry"""
        u_component: MenuSocket
        """Store either the 'Length' or the 'Factor' of the curve as the U component."""
        nucleic_backbone_shape: MenuSocket
        """Nucleic Backbone Shape"""
        nucleic_backbone_radius: FloatSocket
        """Nucleic Backbone Radius"""
        nucleic_backbone_width: FloatSocket
        """Nucleic Backbone Width"""
        nucleic_backbone_thickness: FloatSocket
        """Nucleic Backbone Thickness"""
        base_geometry: GeometrySocket
        """Base Geometry"""
        base_scale: VectorSocket
        """Base Scale"""
        base_resolution: IntegerSocket
        """Base Resolution"""
        base_realize: BooleanSocket
        """Base Realize"""
        shade_smooth: BooleanSocket
        """Apply smooth shading to the created geometry"""
        material: MaterialSocket
        """Material to apply to the resulting geometry"""

    class _Outputs(SocketAccessor):
        geometry: GeometrySocket
        """The generated geometry for the style node group"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        atoms: InputGeometry = None,
        selection: InputBoolean = True,
        quality: InputInteger = 3,
        peptide_radius: InputFloat = 1.6,
        backbone_smoothing: InputFloat = 0.5,
        backbone_threshold: InputFloat = 4.5,
        uv_map: InputBoolean = False,
        u_component: InputMenu | Literal["Factor", "Length"] = "Factor",
        nucleic_backbone_shape: InputMenu
        | Literal["Cicular", "Rectangular"] = "Cicular",
        nucleic_backbone_radius: InputFloat = 2.0,
        nucleic_backbone_width: InputFloat = 4.0,
        nucleic_backbone_thickness: InputFloat = 1.0,
        base_geometry: InputGeometry = None,
        base_scale: InputVector = None,
        base_resolution: InputInteger = 4,
        base_realize: InputBoolean = False,
        shade_smooth: InputBoolean = True,
        material: InputMaterial = None,
    ):
        super().__init__(
            **{
                "Atoms": atoms,
                "Selection": selection,
                "Quality": quality,
                "Peptide Radius": peptide_radius,
                "Backbone Smoothing": backbone_smoothing,
                "Backbone Threshold": backbone_threshold,
                "UV Map": uv_map,
                "U Component": u_component,
                "Nucleic Backbone Shape": nucleic_backbone_shape,
                "Nucleic Backbone Radius": nucleic_backbone_radius,
                "Nucleic Backbone Width": nucleic_backbone_width,
                "Nucleic Backbone Thickness": nucleic_backbone_thickness,
                "Base Geometry": base_geometry,
                "Base Scale": base_scale,
                "Base Resolution": base_resolution,
                "Base Realize": base_realize,
                "Shade Smooth": shade_smooth,
                "Material": material,
            }
        )

    def _build_group(self, tree):
        atoms = tree.inputs.geometry(
            "Atoms", description="Atomic geometry that contains vertices and edges"
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Selection of atoms to apply this style to, discarding unselected points",
            hide_value=True,
        )
        quality = tree.inputs.integer(
            "Quality",
            3,
            description="A lower value results in less geometry, with a higher value meaning better looking but more dense geometry",
            min_value=0,
            max_value=6,
        )
        with tree.inputs.panel("Peptide", default_closed=True):
            peptide_radius = tree.inputs.float("Peptide Radius", 1.6, min_value=0.0)
            with tree.inputs.panel("Backbone"):
                backbone_smoothing = tree.inputs.float(
                    "Backbone Smoothing",
                    0.5,
                    description="Smoothen the sheet ribbons such as beta-sheets",
                    min_value=0.0,
                    max_value=1.0,
                    subtype="FACTOR",
                )
                backbone_threshold = tree.inputs.float(
                    "Backbone Threshold",
                    4.5,
                    description="Distance (Angstroms) over which subsequent CA points are treated as a new chain",
                    min_value=0.0,
                    max_value=10_000.0,
                )
            with tree.inputs.panel(
                "UV Map",
                description="Create a UV map for the generated mesh ribbon.",
                default_closed=True,
            ):
                uv_map = tree.inputs.boolean(
                    "UV Map",
                    False,
                    description="Compute and store the `uv_map` for the final protein ribbon geometry",
                    is_panel_toggle=True,
                )
                u_component = tree.inputs.menu(
                    "U Component",
                    description="Store either the 'Length' or the 'Factor' of the curve as the U component.",
                    expanded=True,
                    optional_label=True,
                )
        with tree.inputs.panel("Nucleic", default_closed=True):
            with tree.inputs.panel("Nucleic Backbone"):
                nucleic_backbone_shape = tree.inputs.menu(
                    "Nucleic Backbone Shape", expanded=True, optional_label=True
                )
                nucleic_backbone_radius = tree.inputs.float(
                    "Nucleic Backbone Radius", 2.0, min_value=0.0
                )
                nucleic_backbone_width = tree.inputs.float(
                    "Nucleic Backbone Width", 4.0, min_value=0.0, max_value=10_000.0
                )
                nucleic_backbone_thickness = tree.inputs.float(
                    "Nucleic Backbone Thickness", 1.0, min_value=0.0, max_value=10_000.0
                )
            with tree.inputs.panel("Nucleic Base", default_closed=True):
                base_geometry = tree.inputs.geometry("Base Geometry")
                base_scale = tree.inputs.vector(
                    "Base Scale",
                    (2.5, 0.5, 7.0),
                    min_value=-10_000.0,
                    max_value=10_000.0,
                )
                base_resolution = tree.inputs.integer(
                    "Base Resolution", 4, min_value=3, max_value=512
                )
                _base_realize = tree.inputs.boolean("Base Realize", False)
        with tree.inputs.panel("Material", default_closed=True):
            shade_smooth = tree.inputs.boolean(
                "Shade Smooth",
                True,
                description="Apply smooth shading to the created geometry",
            )
            material = tree.inputs.material(
                "Material",
                description="Material to apply to the resulting geometry",
                optional_label=True,
            )
        geometry = tree.outputs.geometry(
            "Geometry", description="The generated geometry for the style node group"
        )

        closure_zone = g.ClosureZone()
        atoms_1 = closure_zone.inputs.geometry("Atoms")
        geometry_1 = closure_zone.outputs.geometry("Geometry")
        group = SeparatePolymers(atoms=CheckGeometry(geometry=atoms_1))
        group_1 = MN_utils_style_ribbon_peptide(
            Atoms=group.o.peptide,
            Selection=selection,
            Quality=quality,
            **{"BS Smoothing": backbone_smoothing},
            Radius=peptide_radius,
            **{"Shade Smooth": shade_smooth},
            Material=material,
            **{"UV Map": uv_map, "U Component": u_component},
            Threshold=backbone_threshold,
        )
        menu_switch = g.MenuSwitch.integer(
            nucleic_backbone_shape, {"Cicular": 0, "Rectangular": 1}
        )
        index_switch = g.IndexSwitch.vector(
            menu_switch.o.output,
            (
                (0.0, 0.0, 0.0),
                g.CombineXYZ(
                    y=nucleic_backbone_width, z=nucleic_backbone_thickness, x=1.0
                ),
            ),
        )
        group_2 = MN_utils_style_ribbon_nucleic(
            atoms=group.o.nucleic,
            selection=selection,
            material=material,
            switch=g.IndexSwitch.boolean(menu_switch.o.output, (True, False)),
            backbone_subdivisions=quality * 2,
            backbone_resolution=quality * 4,
            backbone_radius=g.IndexSwitch.float(
                menu_switch.o.output, (nucleic_backbone_radius, 0.0)
            ),
            backbone_shade_smooth=shade_smooth,
            backbone_scale=index_switch,
            base_geometry=base_geometry,
            base_scale=base_scale,
            base_resolution=base_resolution,
        )
        g.JoinGeometry(geometry=(group_1.o.geometry, group_2.o.geometry)) >> geometry_1
        EvaluateOnAtoms(geometry=atoms, closure=closure_zone.closure) >> geometry
        _join_geometry = g.JoinGeometry(geometry=(group_1.o.curve, group_2.o.curve))

        u_component.default_value = "Factor"
        nucleic_backbone_shape.default_value = "Cicular"


ASSET = StyleRibbon

ASSET_METADATA = {
    "catalog_id": "541e6649-2ea6-4225-b1ee-5c0da6f5f1f6",
}
