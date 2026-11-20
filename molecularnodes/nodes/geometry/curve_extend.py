# Node-group asset "Curve Extend" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
import math
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    GeometrySocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputFloat, InputGeometry
from ._shared.control_points_to_curve import ControlPointsToCurve
from ._shared.curve_control_points import CurveControlPoints


class CurveExtend(AssetGeometryGroup):
    """
    Extend the ends of open curves with circular arcs that leave each end with its tangent and curvature, so the curve carries on bending as it was. Each extension adds 4 points that copy the attributes of the end point. Output curves are Bézier with free handles. Extend idea from github.com/Eterea/etr_curves_toolkit (GPL-3.0-or-later, Cristobal Vila)

    Parameters
    ----------
    curve : InputGeometry
        Curve
    selection : InputBoolean
        End points to extend. Cyclic curves are never extended
    start : InputFloat
        Length of the extension before the first point
    end : InputFloat
        Length of the extension after the last point
    follow_curvature : InputBoolean
        Continue the bend of the curve at each end. Off extends straight along the tangent

    Inputs
    ------
    i.curve : GeometrySocket
        Curve
    i.selection : BooleanSocket
        End points to extend. Cyclic curves are never extended
    i.start : FloatSocket
        Length of the extension before the first point
    i.end : FloatSocket
        Length of the extension after the last point
    i.follow_curvature : BooleanSocket
        Continue the bend of the curve at each end. Off extends straight along the tangent

    Outputs
    -------
    o.curve : GeometrySocket
        Curve
    o.extension : BooleanSocket
        Points added by the extension
    """

    _name = "Curve Extend"
    _asset_name = "Curve Extend"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "description": "Extend the ends of open curves with circular arcs that leave each end with its tangent and curvature, so the curve carries on bending as it was. Each extension adds 4 points that copy the attributes of the end point. Output curves are Bézier with free handles. Extend idea from github.com/Eterea/etr_curves_toolkit (GPL-3.0-or-later, Cristobal Vila)"
    }

    class _Inputs(SocketAccessor):
        curve: GeometrySocket
        """Curve"""
        selection: BooleanSocket
        """End points to extend. Cyclic curves are never extended"""
        start: FloatSocket
        """Length of the extension before the first point"""
        end: FloatSocket
        """Length of the extension after the last point"""
        follow_curvature: BooleanSocket
        """Continue the bend of the curve at each end. Off extends straight along the tangent"""

    class _Outputs(SocketAccessor):
        curve: GeometrySocket
        """Curve"""
        extension: BooleanSocket
        """Points added by the extension"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        curve: InputGeometry = None,
        selection: InputBoolean = True,
        start: InputFloat = 0.5,
        end: InputFloat = 0.5,
        follow_curvature: InputBoolean = True,
    ):
        super().__init__(
            **{
                "Curve": curve,
                "Selection": selection,
                "Start": start,
                "End": end,
                "Follow Curvature": follow_curvature,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        curve = tree.inputs.geometry("Curve")
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="End points to extend. Cyclic curves are never extended",
            hide_value=True,
        )
        start = tree.inputs.float(
            "Start",
            0.5,
            description="Length of the extension before the first point",
            min_value=0.0,
            subtype="DISTANCE",
        )
        end = tree.inputs.float(
            "End",
            0.5,
            description="Length of the extension after the last point",
            min_value=0.0,
            subtype="DISTANCE",
        )
        follow_curvature = tree.inputs.boolean(
            "Follow Curvature",
            True,
            description="Continue the bend of the curve at each end. Off extends straight along the tangent",
        )
        curve_1 = tree.outputs.geometry("Curve")
        extension = tree.outputs.boolean(
            "Extension", description="Points added by the extension"
        )

        set_spline_type = g.SetSplineType.bezier(curve)
        with g.Frame("End frame"):
            curve_handle_positions = g.CurveHandlePositions()
            curve_handle_positions_1 = g.CurveHandlePositions()
            curve_tangent = g.CurveTangent()
            curve_of_point = g.CurveOfPoint()
            points_of_curve = g.PointsOfCurve(curve_of_point.o.curve_index)
            compare = g.Compare.integer.equal(
                curve_of_point.o.index_in_curve, points_of_curve.o.total - 1
            )
            boolean_math = (
                selection
                & ~g.IsSplineCyclic().o.cyclic
                & (points_of_curve.o.total >= 2)
                & (
                    g.Compare.integer.equal(curve_of_point.o.index_in_curve, 0).o.result
                    | compare
                )
            )
            position = g.Position()
            _string = g.String(
                "The last segment at each end, read outward, has first and second derivatives d1 = 3 (P3 - P2) and d2 = 6 (P3 - 2 P2 + P1). The curvature vector is the part of d2 across the tangent over |d1|^2. The arc is capped at one full turn and split into 4 pieces of at most 90 degrees, each a cubic with handles 4/3 tan(angle / 4) / curvature long."
            )
            switch = compare.o.result.switch.integer(
                g.Index().o.index + 1, g.Index().o.index - 1
            )
            switch_1 = compare.o.result.switch.vector(
                curve_handle_positions_1.o.right, curve_handle_positions.o.left
            )
            switch_2 = compare.o.result.switch.vector(
                curve_handle_positions.o.left.point.at(switch),
                curve_handle_positions_1.o.right.point.at(switch),
            )
            switch_3 = compare.o.result.switch.float(start, end)
            vector_math = (position.o.position - switch_1) * 3.0
            vector_math_1 = (position.o.position - switch_1 * 2.0 + switch_2) * 6.0
            switch_4 = compare.o.result.switch.vector(
                curve_tangent.o.tangent * -1.0, curve_tangent
            )
            vector_math_2 = vector_math.dot(vector_math)
            vector_math_3 = (
                vector_math_1 - switch_4 * vector_math_1.dot(switch_4)
            ) / vector_math_2.max(0.000000000001)
            switch_5 = (
                follow_curvature & (vector_math_2 > 0.000000000001)
            ).switch.vector((0.0, 0.0, 0.0), vector_math_3)
            vector_math_4 = switch_5.length()
            compare_1 = vector_math_4 > 0.000001
            math_1 = (
                compare_1.switch.float(switch_3, switch_3.min(math.tau / vector_math_4))
                / 4.0
            )
            switch_6 = compare_1.switch.float(
                math_1 / 3.0,
                (vector_math_4 * math_1 / 4.0).tan() * 4.0 / 3.0 / vector_math_4,
            )
            capture = g.CaptureAttribute.point(geometry=set_spline_type)
            tangent = capture.items.vector(switch_4, "Tangent")
            normal = capture.items.vector(switch_5.normalize(), "Normal")
            curvature = capture.items.float(vector_math_4, "Curvature")
            step = capture.items.float(math_1, "Step")
            handle = capture.items.float(switch_6, "Handle")
            extend = capture.items.boolean(boolean_math & (switch_3 > 0.0), "Extend")
            is_end = capture.items.boolean(compare, "Is End")
        vector_math_5 = position.o.position + tangent.output * handle.output
        set_handle_positions = (
            g.SetHandleType.free(capture.o.geometry, extend.output)
            >> g.SetHandlePositions(
                selection=extend.output & ~is_end.output, position=vector_math_5
            )
            >> g.SetHandlePositions.right(
                selection=extend.output & is_end.output, position=vector_math_5
            )
        )
        curve_control_points = CurveControlPoints(curve=set_handle_positions)
        capture_1 = g.CaptureAttribute.point(geometry=curve_control_points)
        order = capture_1.items.float(g.Index(), "Order")
        with g.Frame("Arc points"):
            _string_1 = g.String(
                "Each extended end is copied 4 times. Copy j sits at arc length s = (j + 1) step from the end: P + t sin(k s) / k + n (1 - cos(k s)) / k, with tangent t cos(k s) + n sin(k s). Copies before the first point sort before it."
            )
            compare_2 = curvature.output > 0.000001
            switch_7 = is_end.output.switch.float(-1.0, 1.0)
            duplicate_elements = g.DuplicateElements.point(
                capture_1.o.geometry, extend.output, 4
            )
            integer_math = duplicate_elements.o.duplicate_index + 1
            math_2 = integer_math * step.output
            math_3 = curvature.output * math_2
            vector_math_6 = (
                position.o.position
                + tangent.output
                * compare_2.switch.float(math_2, math_3.sin() / curvature.output)
            )
            vector_math_7 = vector_math_6 + normal.output * compare_2.switch.float(
                true=(1.0 - math_3.cos()) / curvature.output
            )
            vector_math_8 = (
                tangent.output * math_3.cos() + normal.output * math_3.sin()
            ) * (handle.output * switch_7)
            capture_2 = g.CaptureAttribute.point(geometry=duplicate_elements)
            position_1 = capture_2.items.vector(vector_math_7, "Position")
            handle_left = capture_2.items.vector(
                vector_math_7 - vector_math_8, "Handle Left"
            )
            handle_right = capture_2.items.vector(
                vector_math_7 + vector_math_8, "Handle Right"
            )
            offset = capture_2.items.float(integer_math * switch_7 / 5.0, "Offset")
            extension_1 = capture_2.items.boolean(True, "Extension")
        join_geometry = g.JoinGeometry(
            (
                capture_1.o.geometry,
                capture_2.o.geometry >> g.SetPosition(position=position_1.output),
            )
        )
        switch_8 = extension_1.output.switch.vector(
            curve_control_points.o.handle_right, handle_right.output
        )
        (
            ControlPointsToCurve(
                points=join_geometry,
                curve_group_id=curve_control_points.o.curve_index,
                weight=order.output + offset.output,
                handle_left=extension_1.output.switch.vector(
                    curve_control_points.o.handle_left, handle_left.output
                ),
                handle_right=switch_8,
                cyclic=curve_control_points.o.cyclic,
            )
            >> curve_1
        )

        extension_1.output >> extension


ASSET = CurveExtend

ASSET_METADATA = {
    "description": "Extend the ends of open curves with circular arcs that leave each end with its tangent and curvature, so the curve carries on bending as it was. Each extension adds 4 points that copy the attributes of the end point. Output curves are Bézier with free handles. Extend idea from github.com/Eterea/etr_curves_toolkit (GPL-3.0-or-later, Cristobal Vila)",
    "catalog_id": "9c167a5c-d0a6-457d-9e9c-90f3edd29e10",
}
