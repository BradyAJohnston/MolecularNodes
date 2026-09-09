# Node group 'Curve Split Splines' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from nodebpy import geometry as g
from nodebpy.builder import (
    BooleanSocket,
    CustomGeometryGroup,
    FloatSocket,
    GeometrySocket,
    IntegerSocket,
    MenuSocket,
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
)
from ..angstrom_to_world import AngstromToWorld
from ..offset_point_along_curve import OffsetPointAlongCurve
from ..offset_vector import OffsetVector


class CurveSplitSplines(CustomGeometryGroup):
    """
    Spline the given curves into new separate splines based on their Curve Group ID, additionally dropping any non-selected points

    Parameters
    ----------
    curve : InputGeometry
        Curve
    selection : InputBoolean
        Selection
    curve_normal : InputMenu | Literal["Minimum Twist", "Free"]
        Curve Normal
    curve_group_id : InputInteger
        Curve Group ID
    distance_split : InputMenu | Literal["Ignore Distance", "Split Distance"]
        Distance Split
    distance_cutoff : InputFloat
        Distance Cutoff
    rotation : InputRotation
        Rotation
    offset_amount : InputFloat
        Offset Amount
    offset_spline_type : InputMenu | Literal["Poly", "Bezier"]
        Offset Spline Type
    offset_resolution : InputInteger
        Offset Resolution

    Inputs
    ------
    i.curve : GeometrySocket
        Curve
    i.selection : BooleanSocket
        Selection
    i.curve_normal : MenuSocket
        Curve Normal
    i.curve_group_id : IntegerSocket
        Curve Group ID
    i.distance_split : MenuSocket
        Distance Split
    i.distance_cutoff : FloatSocket
        Distance Cutoff
    i.rotation : RotationSocket
        Rotation
    i.offset_amount : FloatSocket
        Offset Amount
    i.offset_spline_type : MenuSocket
        Offset Spline Type
    i.offset_resolution : IntegerSocket
        Offset Resolution

    Outputs
    -------
    o.curve : GeometrySocket
        Curve
    o.index : IntegerSocket
        Index
    o.rotation : RotationSocket
        Rotation
    o.offset_rotation : RotationSocket
        Offset Rotation
    o.offset_position : VectorSocket
        Offset Position
    o.offset_tangent : VectorSocket
        Offset Tangent
    o.offset_normal : VectorSocket
        Offset Normal
    """

    _name = "Curve Split Splines"
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "description": "Spline the given curves into new separate splines based on their Curve Group ID, additionally dropping any non-selected points",
        "node_tool_idname": "geometry.curve_split_splines",
    }

    class _Inputs(SocketAccessor):
        curve: GeometrySocket
        """Curve"""
        selection: BooleanSocket
        """Selection"""
        curve_normal: MenuSocket
        """Curve Normal"""
        curve_group_id: IntegerSocket
        """Curve Group ID"""
        distance_split: MenuSocket
        """Distance Split"""
        distance_cutoff: FloatSocket
        """Distance Cutoff"""
        rotation: RotationSocket
        """Rotation"""
        offset_amount: FloatSocket
        """Offset Amount"""
        offset_spline_type: MenuSocket
        """Offset Spline Type"""
        offset_resolution: IntegerSocket
        """Offset Resolution"""

    class _Outputs(SocketAccessor):
        curve: GeometrySocket
        """Curve"""
        index: IntegerSocket
        """Index"""
        rotation: RotationSocket
        """Rotation"""
        offset_rotation: RotationSocket
        """Offset Rotation"""
        offset_position: VectorSocket
        """Offset Position"""
        offset_tangent: VectorSocket
        """Offset Tangent"""
        offset_normal: VectorSocket
        """Offset Normal"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        curve: InputGeometry = None,
        selection: InputBoolean = True,
        curve_normal: InputMenu | Literal["Minimum Twist", "Free"] = "Free",
        curve_group_id: InputInteger = 0,
        distance_split: InputMenu
        | Literal["Ignore Distance", "Split Distance"] = "Split Distance",
        distance_cutoff: InputFloat = 0.0,
        rotation: InputRotation = None,
        offset_amount: InputFloat = 0.0,
        offset_spline_type: InputMenu | Literal["Poly", "Bezier"] = "Bezier",
        offset_resolution: InputInteger = 12,
    ):
        super().__init__(
            **{
                "Curve": curve,
                "Selection": selection,
                "Curve Normal": curve_normal,
                "Curve Group ID": curve_group_id,
                "Distance Split": distance_split,
                "Distance Cutoff": distance_cutoff,
                "Rotation": rotation,
                "Offset Amount": offset_amount,
                "Offset Spline Type": offset_spline_type,
                "Offset Resolution": offset_resolution,
            }
        )

    def _build_group(self, tree):
        curve = tree.inputs.geometry("Curve")
        selection = tree.inputs.boolean("Selection", True, hide_value=True)
        curve_normal = tree.inputs.menu("Curve Normal", optional_label=True)
        curve_group_id = tree.inputs.integer("Curve Group ID", 0)
        distance_split = tree.inputs.menu("Distance Split", optional_label=True)
        distance_cutoff = tree.inputs.float(
            "Distance Cutoff", 0.0, min_value=-10_000.0, max_value=10_000.0
        )
        rotation = tree.inputs.rotation("Rotation", (0.0, 0.0, 0.0), hide_value=True)
        with tree.panel("Offset"):
            offset_amount = tree.inputs.float(
                "Offset Amount", 0.0, min_value=-10_000.0, max_value=10_000.0
            )
            offset_spline_type = tree.inputs.menu(
                "Offset Spline Type", optional_label=True
            )
            offset_resolution = tree.inputs.integer(
                "Offset Resolution", 12, min_value=1
            )
        curve_1 = tree.outputs.geometry("Curve")
        index = tree.outputs.integer("Index")
        rotation_1 = tree.outputs.rotation("Rotation")
        with tree.panel("Offset"):
            offset_rotation = tree.outputs.rotation("Offset Rotation")
            offset_position = tree.outputs.vector("Offset Position")
            offset_tangent = tree.outputs.vector("Offset Tangent")
            offset_normal = tree.outputs.vector("Offset Normal")

        curve_handle_positions = g.CurveHandlePositions()
        capture = g.CaptureAttribute.point(geometry=curve)
        factor = capture.items.float(
            "Factor", OffsetPointAlongCurve(offset=offset_amount).o.factor
        )
        curve_index = capture.items.integer(
            "Curve Index", g.CurveOfPoint().o.curve_index
        )
        selection_1 = capture.items.boolean("Selection", selection)
        capture.items.integer("Curve Group ID", curve_group_id)
        distance_cutoff_1 = capture.items.float("Distance Cutoff", distance_cutoff)
        rotation_2 = capture.items.rotation("Rotation", rotation)
        index_1 = capture.items.integer("Index", g.Index())
        normal = capture.items.vector(
            "Normal", g.Normal(legacy_corner_normals=True).o.normal
        )
        position = capture.items.vector("Position", g.Position())
        trailing = capture.items.integer(
            "Trailing", g.AccumulateField.point.integer(~selection).o.trailing
        )
        with g.Frame("Potentially sample a new interpolated point along the curve"):
            sample_curve = g.MenuSwitch.geometry(
                offset_spline_type,
                {
                    "Poly": capture.o.geometry,
                    "Bezier": g.SetSplineResolution(
                        geometry=capture.o.geometry, resolution=offset_resolution
                    ),
                },
            ) >> g.SampleCurve(
                value=rotation_2.output,
                factor=factor.output,
                curve_index=curve_index.output,
                data_type="QUATERNION",
            )
        set_spline_resolution = g.SetHandleType(
            curve=g.SetSplineType.bezier(capture.o.geometry)
        ) >> g.SetSplineResolution(resolution=1)
        capture_1 = g.CaptureAttribute.point(geometry=set_spline_resolution)
        left = capture_1.items.vector("Left", curve_handle_positions.o.left)
        right = capture_1.items.vector("Right", curve_handle_positions.o.right)
        separate_geometry = (
            capture_1.o.geometry
            >> g.CurveToPoints.evaluated()
            >> g.SeparateGeometry.point(selection=selection_1.output)
        )
        compare = (
            AngstromToWorld(
                angstrom=OffsetVector(offset=-1).o.value.distance(g.Position())
            )
            > distance_cutoff_1.output
        )
        math_1 = g.Math(
            value=trailing.output,
            value_001=g.AccumulateField.point.integer(compare).o.leading,
        )
        menu_switch = g.MenuSwitch.integer(
            distance_split,
            {"Ignore Distance": trailing.output, "Split Distance": math_1},
        )
        points_to_curves = g.PointsToCurves(
            points=g.SetPosition(
                geometry=separate_geometry.o.selection, position=position.output
            ),
            curve_group_id=menu_switch.o.output,
        )
        set_spline_type = (
            g.MenuSwitch.geometry(
                curve_normal,
                {
                    "Minimum Twist": points_to_curves,
                    "Free": g.SetCurveNormal(
                        curve=points_to_curves, normal=normal.output, mode="Free"
                    ),
                },
            )
            >> g.SetSplineType.bezier()
        )
        (
            g.SetHandleType(curve=set_spline_type)
            >> g.SetHandlePositions(position=left.output)
            >> g.SetHandlePositions.right(position=right.output)
            >> curve_1
        )
        viewer = g.Viewer()
        separate_geometry >> viewer

        index_1.output >> index
        rotation_2.output >> rotation_1
        sample_curve >> offset_rotation
        sample_curve.o.position >> offset_position
        sample_curve.o.tangent >> offset_tangent
        sample_curve.o.normal >> offset_normal

        curve_normal.default_value = "Free"
        distance_split.default_value = "Split Distance"
        offset_spline_type.default_value = "Bezier"
