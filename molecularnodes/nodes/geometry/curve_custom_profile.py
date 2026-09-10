# Node-group asset "Curve Custom Profile" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
import math
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
    MenuSocket,
    PackageLibrary,
    RotationSocket,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import (
    InputBoolean,
    InputFloat,
    InputGeometry,
    InputInteger,
    InputMenu,
    InputRotation,
    InputVector,
)
from ._shared.check_end_face_corner import CheckEndFaceCorner
from ._shared.mn_units import MNUnits
from ._shared.profile_type_picker import ProfileTypePicker
from .curve_rotation import CurveRotation
from .fallback_geometry import FallbackGeometry
from .sample_position import SamplePosition


class CurveCustomProfile(AssetGeometryGroup):
    """
    Curve Custom Profile

    Parameters
    ----------
    curve : InputGeometry
        Curve
    subdivisions : InputInteger
        Subdivisions
    profile_type : InputMenu | Literal["Default Profile", "Custom Profile"]
        Profile Type
    uv_map : InputBoolean
        Compute and store the `uv_map` attribute on the `Face Corner` domain of the final mesh
    u_component : InputMenu | Literal["Factor", "Length"]
        U Component
    socket_6 : InputRotation
        Profile Rotation
    profile_scale : InputVector
        Profile Scale
    profile_curve : InputGeometry
        Profile Curve
    profile_resolution : InputInteger
        Profile Resolution
    profile_radius : InputFloat
        Profile Radius
    input_14 : InputFloat
        Profile Rotation

    Inputs
    ------
    i.curve : GeometrySocket
        Curve
    i.subdivisions : IntegerSocket
        Subdivisions
    i.profile_type : MenuSocket
        Profile Type
    i.uv_map : BooleanSocket
        Compute and store the `uv_map` attribute on the `Face Corner` domain of the final mesh
    i.u_component : MenuSocket
        U Component
    i.socket_6 : RotationSocket
        Profile Rotation
    i.profile_scale : VectorSocket
        Profile Scale
    i.profile_curve : GeometrySocket
        Profile Curve
    i.profile_resolution : IntegerSocket
        Profile Resolution
    i.profile_radius : FloatSocket
        Profile Radius
    i.input_14 : FloatSocket
        Profile Rotation

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "Curve Custom Profile"
    _asset_name = "Curve Custom Profile"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry.curve_custom_profile"}

    class _Inputs(SocketAccessor):
        curve: GeometrySocket
        """Curve"""
        subdivisions: IntegerSocket
        """Subdivisions"""
        profile_type: MenuSocket
        """Profile Type"""
        uv_map: BooleanSocket
        """Compute and store the `uv_map` attribute on the `Face Corner` domain of the final mesh"""
        u_component: MenuSocket
        """U Component"""
        socket_6: RotationSocket
        """Profile Rotation"""
        profile_scale: VectorSocket
        """Profile Scale"""
        profile_curve: GeometrySocket
        """Profile Curve"""
        profile_resolution: IntegerSocket
        """Profile Resolution"""
        profile_radius: FloatSocket
        """Profile Radius"""
        input_14: FloatSocket
        """Profile Rotation"""

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
        curve: InputGeometry = None,
        subdivisions: InputInteger = 6,
        profile_type: InputMenu
        | Literal["Default Profile", "Custom Profile"] = "Custom Profile",
        uv_map: InputBoolean = False,
        u_component: InputMenu | Literal["Factor", "Length"] = "Factor",
        socket_6: InputRotation = None,
        profile_scale: InputVector = None,
        profile_curve: InputGeometry = None,
        profile_resolution: InputInteger = 4,
        profile_radius: InputFloat = 1.0,
        input_14: InputFloat = math.pi / 4,
    ):
        super().__init__(
            **{
                "Curve": curve,
                "Subdivisions": subdivisions,
                "Profile Type": profile_type,
                "UV Map": uv_map,
                "U Component": u_component,
                "Profile Scale": profile_scale,
                "Profile Curve": profile_curve,
                "Profile Resolution": profile_resolution,
                "Profile Radius": profile_radius,
            },
            _named_links=[
                ("Profile Rotation", socket_6),
                ("Profile Rotation", input_14),
            ],
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        curve = tree.inputs.geometry("Curve")
        subdivisions = tree.inputs.integer("Subdivisions", 6, min_value=1)
        profile_type = tree.inputs.menu("Profile Type", optional_label=True)
        uv_map = tree.inputs.boolean(
            "UV Map",
            False,
            description="Compute and store the `uv_map` attribute on the `Face Corner` domain of the final mesh",
        )
        u_component = tree.inputs.menu("U Component", optional_label=True)
        with tree.inputs.panel("Profile", default_closed=True):
            profile_rotation = tree.inputs.rotation("Profile Rotation", (0.0, 0.0, 0.0))
            profile_scale = tree.inputs.vector(
                "Profile Scale", (1.0, 1.0, 1.0), subtype="XYZ"
            )
            profile_curve = tree.inputs.geometry("Profile Curve")
            profile_resolution = tree.inputs.integer(
                "Profile Resolution", 4, min_value=3, max_value=512
            )
            profile_radius = tree.inputs.float(
                "Profile Radius", 1.0, min_value=0.0, subtype="DISTANCE"
            )
            profile_rotation_1 = tree.inputs.float(
                "Profile Rotation", math.pi / 4, min_value=-10_000.0, max_value=10_000.0
            )
        geometry = tree.outputs.geometry("Geometry")

        curve_circle = g.CurveCircle(
            resolution=profile_resolution,
            radius=MNUnits(value=profile_radius).o.angstrom,
        )
        index_switch = g.IndexSwitch.rotation(
            ProfileTypePicker(menu=profile_type), (CurveRotation(), profile_rotation)
        )
        spline_parameter = g.SplineParameter()
        capture = g.CaptureAttribute.point(geometry=curve)
        rotation = capture.items.rotation("Rotation", index_switch)
        scale = capture.items.vector("Scale", profile_scale)
        factor = capture.items.float("Factor", spline_parameter.o.factor)
        length = capture.items.float("Length", spline_parameter.o.length)
        index = capture.items.integer("Index", spline_parameter.o.index)
        resample_curve = (
            capture.o.geometry
            >> g.SetSplineResolution(resolution=subdivisions)
            >> g.ResampleCurve(mode="Evaluated", length=0.1)
        )
        spline_parameter_1 = g.SplineParameter()
        switch = CheckEndFaceCorner(
            captured_index=index.output
        ).o.is_end_face_corner.switch.float(
            g.MenuSwitch.float(
                u_component, {"Factor": factor.output, "Length": length.output}
            ).o.output,
            g.IndexSwitch.float(items=(1.0, 38.0)),
        )
        transform_geometry = g.TransformGeometry(
            geometry=curve_circle,
            rotation=g.AxisAngleToRotation(
                angle=profile_rotation_1 + g.Math.to_radians(45.0)
            ),
        )
        capture_1 = g.CaptureAttribute.point(
            geometry=FallbackGeometry(
                geometry=profile_curve, fallback=transform_geometry
            )
        )
        factor_1 = capture_1.items.float("Factor", spline_parameter_1.o.factor)
        index_1 = capture_1.items.integer("Index", spline_parameter_1.o.index)
        curve_to_mesh = g.CurveToMesh(
            curve=resample_curve,
            profile_curve=capture_1.o.geometry,
            scale=g.Radius(),
            fill_caps=True,
        )
        instance_on_points = g.InstanceOnPoints(
            points=resample_curve,
            instance=capture_1.o.geometry,
            rotation=rotation.output,
            scale=scale.output,
        )
        group = SamplePosition(
            geometry=g.RealizeInstances(
                geometry=instance_on_points, realize_to_point_domain=True
            )
        )
        switch_1 = CheckEndFaceCorner(
            captured_index=index_1.output
        ).o.is_end_face_corner.switch.float(factor_1.output, 1.0)
        store_named_attribute = g.StoreNamedAttribute.corner.vector_2d(
            curve_to_mesh, name="uv_map", value=g.CombineXYZ(x=switch, y=switch_1)
        )
        switch_2 = uv_map.switch.geometry(curve_to_mesh, store_named_attribute)
        (
            g.IndexSwitch.geometry(
                ProfileTypePicker(menu=profile_type),
                (switch_2, g.SetPosition(geometry=switch_2, position=group)),
            )
            >> geometry
        )

        profile_type.default_value = "Custom Profile"
        u_component.default_value = "Factor"


ASSET = CurveCustomProfile

ASSET_METADATA = {
    "catalog_id": "9c167a5c-d0a6-457d-9e9c-90f3edd29e10",
}
