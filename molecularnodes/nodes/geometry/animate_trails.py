# Node-group asset "Animate Trails" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    GeometrySocket,
    IntegerSocket,
    MaterialSocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import (
    InputBoolean,
    InputFloat,
    InputGeometry,
    InputInteger,
    InputMaterial,
    InputMenu,
)
from .break_curves import BreakCurves
from .lag_geometry import LagGeometry
from .vdw_radii import VDWRadii


class AnimateTrails(AssetGeometryGroup):
    """
    Animate Trails

    Parameters
    ----------
    atoms : InputGeometry
        Atomic geometry that contains vertices and edges
    selection : InputBoolean
        Selection of atoms to apply this node to
    trail_threshold : InputMenu | Literal["Unlimited", "Threshold"]
        Trail Threshold
    trail_curve_type : InputMenu | Literal["Poly", "Bezier"]
        Trail Curve Type
    trail_frames : InputInteger
        Number of previous frames from the trajectory to display
    trail_radius : InputFloat
        Trail Radius
    trail_resolution : InputInteger
        Tail radial resolution
    trail_subdivisions : InputInteger
        Trail Subdivisions
    trail_cutoff : InputFloat
        Threshold over which the spline is broken up
    shade_smooth : InputBoolean
        Apply smooth shading to the created geometry
    material : InputMaterial
        Material to apply to the resulting geometry

    Inputs
    ------
    i.atoms : GeometrySocket
        Atomic geometry that contains vertices and edges
    i.selection : BooleanSocket
        Selection of atoms to apply this node to
    i.trail_threshold : MenuSocket
        Trail Threshold
    i.trail_curve_type : MenuSocket
        Trail Curve Type
    i.trail_frames : IntegerSocket
        Number of previous frames from the trajectory to display
    i.trail_radius : FloatSocket
        Trail Radius
    i.trail_resolution : IntegerSocket
        Tail radial resolution
    i.trail_subdivisions : IntegerSocket
        Trail Subdivisions
    i.trail_cutoff : FloatSocket
        Threshold over which the spline is broken up
    i.shade_smooth : BooleanSocket
        Apply smooth shading to the created geometry
    i.material : MaterialSocket
        Material to apply to the resulting geometry

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "Animate Trails"
    _asset_name = "Animate Trails"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry.animate_trails"}

    class _Inputs(SocketAccessor):
        atoms: GeometrySocket
        """Atomic geometry that contains vertices and edges"""
        selection: BooleanSocket
        """Selection of atoms to apply this node to"""
        trail_threshold: MenuSocket
        """Trail Threshold"""
        trail_curve_type: MenuSocket
        """Trail Curve Type"""
        trail_frames: IntegerSocket
        """Number of previous frames from the trajectory to display"""
        trail_radius: FloatSocket
        """Trail Radius"""
        trail_resolution: IntegerSocket
        """Tail radial resolution"""
        trail_subdivisions: IntegerSocket
        """Trail Subdivisions"""
        trail_cutoff: FloatSocket
        """Threshold over which the spline is broken up"""
        shade_smooth: BooleanSocket
        """Apply smooth shading to the created geometry"""
        material: MaterialSocket
        """Material to apply to the resulting geometry"""

    class _Outputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        atoms: InputGeometry = None,
        selection: InputBoolean = True,
        trail_threshold: InputMenu | Literal["Unlimited", "Threshold"] = "Unlimited",
        trail_curve_type: InputMenu | Literal["Poly", "Bezier"] = "Poly",
        trail_frames: InputInteger = 5,
        trail_radius: InputFloat = 1.0,
        trail_resolution: InputInteger = 6,
        trail_subdivisions: InputInteger = 6,
        trail_cutoff: InputFloat = 10.0,
        shade_smooth: InputBoolean = True,
        material: InputMaterial = None,
    ):
        super().__init__(
            **{
                "Atoms": atoms,
                "Selection": selection,
                "Trail Threshold": trail_threshold,
                "Trail Curve Type": trail_curve_type,
                "Trail Frames": trail_frames,
                "Trail Radius": trail_radius,
                "Trail Resolution": trail_resolution,
                "Trail Subdivisions": trail_subdivisions,
                "Trail Cutoff": trail_cutoff,
                "Shade Smooth": shade_smooth,
                "Material": material,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atoms = tree.inputs.geometry(
            "Atoms", description="Atomic geometry that contains vertices and edges"
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Selection of atoms to apply this node to",
            hide_value=True,
        )
        with tree.inputs.panel("Trail"):
            trail_threshold = tree.inputs.menu(
                "Trail Threshold", expanded=True, optional_label=True
            )
            trail_curve_type = tree.inputs.menu(
                "Trail Curve Type", expanded=True, optional_label=True
            )
            trail_frames = tree.inputs.integer(
                "Trail Frames",
                5,
                description="Number of previous frames from the trajectory to display",
                min_value=1,
                max_value=100000,
            )
            trail_radius = tree.inputs.float(
                "Trail Radius", 1.0, min_value=0.0, max_value=10_000.0
            )
            trail_resolution = tree.inputs.integer(
                "Trail Resolution",
                6,
                description="Tail radial resolution",
                min_value=3,
                max_value=32,
            )
            trail_subdivisions = tree.inputs.integer(
                "Trail Subdivisions", 6, min_value=1, max_value=16
            )
            trail_cutoff = tree.inputs.float(
                "Trail Cutoff",
                10.0,
                description="Threshold over which the spline is broken up",
                min_value=0.0,
                max_value=10_000.0,
                subtype="DISTANCE",
            )
        with tree.inputs.panel("Material"):
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
        geometry = tree.outputs.geometry("Geometry")

        group = LagGeometry(input=atoms, selection=selection, count=trail_frames)
        reverse_curve = (
            group
            >> g.MeshToPoints(radius=0.05)
            >> g.PointsToCurves(curve_group_id=group.o.index, weight=group.o.lag_index)
            >> g.ReverseCurve()
        )
        capture = g.CaptureAttribute.point(geometry=reverse_curve)
        factor = capture.items.float("Factor", g.SplineParameter().o.factor)
        menu_switch = g.MenuSwitch.geometry(
            trail_threshold,
            {
                "Unlimited": capture.o.geometry,
                "Threshold": BreakCurves(
                    curves=capture.o.geometry, threshold=trail_cutoff
                ),
            },
        )
        set_spline_resolution = g.SetHandleType(
            curve=g.SetSplineType.bezier(menu_switch)
        ) >> g.SetSplineResolution(resolution=trail_subdivisions)
        (
            g.MenuSwitch.geometry(
                trail_curve_type, {"Poly": menu_switch, "Bezier": set_spline_resolution}
            )
            >> g.CurveToMesh(
                profile_curve=g.CurveCircle(resolution=trail_resolution),
                scale=factor.output * (VDWRadii().o.vdw_radii * trail_radius),
                fill_caps=True,
            )
            >> g.SetShadeSmooth.face(shade_smooth=shade_smooth)
            >> g.SetMaterial(material=material)
            >> geometry
        )

        trail_threshold.default_value = "Unlimited"
        trail_curve_type.default_value = "Poly"


ASSET = AnimateTrails

ASSET_METADATA = {
    "catalog_id": "85730213-4c2e-469f-b333-52ac53adf274",
}

DATABLOCK_DEPENDENCIES = {
    "materials": ("MN Ambient Occlusion",),
}
