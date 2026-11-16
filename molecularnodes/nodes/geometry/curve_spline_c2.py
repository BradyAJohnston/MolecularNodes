# Node-group asset "Curve Spline C2" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    GeometrySocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputGeometry, InputInteger


class CurveSplineC2(AssetGeometryGroup):
    """
    Convert curves to Bézier curves that pass through every control point with continuous curvature (a C2 cubic spline). Open curves have zero curvature at their ends. Curves with fewer than 3 points are converted to Bézier but not solved. Chord Length spaces the spline by the distance between points. Port of SolveCubicC2 from github.com/aurimasg/splines (MIT, Aurimas Gasiulis)

    Parameters
    ----------
    curve : InputGeometry
        Curves whose points the spline passes through
    chord_length : InputBoolean
        Space the spline parameter by the distance between points, which stops overshoot and loops where spacing is uneven. Off spaces every segment equally, as the reference implementation does
    iterations : InputInteger
        Solver steps. Each step at least halves the error, so 24 reaches float precision for any curve length

    Inputs
    ------
    i.curve : GeometrySocket
        Curves whose points the spline passes through
    i.chord_length : BooleanSocket
        Space the spline parameter by the distance between points, which stops overshoot and loops where spacing is uneven. Off spaces every segment equally, as the reference implementation does
    i.iterations : IntegerSocket
        Solver steps. Each step at least halves the error, so 24 reaches float precision for any curve length

    Outputs
    -------
    o.curve : GeometrySocket
        Bézier curves with free handles set to the C2 spline
    """

    _name = "Curve Spline C2"
    _asset_name = "Curve Spline C2"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "description": "Convert curves to Bézier curves that pass through every control point with continuous curvature (a C2 cubic spline). Open curves have zero curvature at their ends. Curves with fewer than 3 points are converted to Bézier but not solved. Chord Length spaces the spline by the distance between points. Port of SolveCubicC2 from github.com/aurimasg/splines (MIT, Aurimas Gasiulis)"
    }

    class _Inputs(SocketAccessor):
        curve: GeometrySocket
        """Curves whose points the spline passes through"""
        chord_length: BooleanSocket
        """Space the spline parameter by the distance between points, which stops overshoot and loops where spacing is uneven. Off spaces every segment equally, as the reference implementation does"""
        iterations: IntegerSocket
        """Solver steps. Each step at least halves the error, so 24 reaches float precision for any curve length"""

    class _Outputs(SocketAccessor):
        curve: GeometrySocket
        """Bézier curves with free handles set to the C2 spline"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        curve: InputGeometry = None,
        chord_length: InputBoolean = True,
        iterations: InputInteger = 24,
    ):
        super().__init__(
            **{"Curve": curve, "Chord Length": chord_length, "Iterations": iterations}
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        curve = tree.inputs.geometry(
            "Curve", description="Curves whose points the spline passes through"
        )
        chord_length = tree.inputs.boolean(
            "Chord Length",
            True,
            description="Space the spline parameter by the distance between points, which stops overshoot and loops where spacing is uneven. Off spaces every segment equally, as the reference implementation does",
        )
        iterations = tree.inputs.integer(
            "Iterations",
            24,
            description="Solver steps. Each step at least halves the error, so 24 reaches float precision for any curve length",
            min_value=0,
            max_value=1000,
        )
        curve_1 = tree.outputs.geometry(
            "Curve", description="Bézier curves with free handles set to the C2 spline"
        )

        curve_of_point = g.CurveOfPoint()
        is_spline_cyclic = g.IsSplineCyclic()
        boolean_math = ~is_spline_cyclic.o.cyclic & g.Compare.integer.equal(
            curve_of_point.o.index_in_curve, 0
        )
        offset_point_in_curve = g.OffsetPointInCurve(offset=-1)
        offset_point_in_curve_1 = g.OffsetPointInCurve(offset=1)
        points_of_curve = g.PointsOfCurve(curve_of_point.o.curve_index)
        compare = points_of_curve.o.total >= 3
        boolean_math_1 = ~is_spline_cyclic.o.cyclic & g.Compare.integer.equal(
            curve_of_point.o.index_in_curve, points_of_curve.o.total - 1
        )
        position = g.Position()
        evaluate_at_index = position.o.position.point.at(
            offset_point_in_curve.o.point_index
        )
        evaluate_at_index_1 = position.o.position.point.at(
            offset_point_in_curve_1.o.point_index
        )
        set_handle_positions = g.SetHandleType.free(
            g.SetSplineType.bezier(curve)
        ) >> g.SetHandlePositions.right(selection=compare, position=position)
        with g.Frame("Parameter spacing"):
            _string = g.String(
                "h_prev and h_next are the parameter lengths of the segments before and after each point: the chord length, or 1 for uniform spacing. Open ends copy the one segment they have, and a floor stops coincident points dividing by zero."
            )
            vector_math = (position.o.position - evaluate_at_index).length()
            vector_math_1 = (evaluate_at_index_1 - position).length()
            math_1 = offset_point_in_curve.o.is_valid_offset.switch.float(
                vector_math_1, vector_math
            ).max(0.000001)
            switch = chord_length.switch.float(1.0, math_1)
            math_2 = offset_point_in_curve_1.o.is_valid_offset.switch.float(
                vector_math, vector_math_1
            ).max(0.000001)
            switch_1 = chord_length.switch.float(1.0, math_2)
        with g.Frame("Solve tangents"):
            repeat_zone = g.RepeatZone(iterations)
            geometry = repeat_zone.items.geometry(set_handle_positions, "Geometry")
            _string_1 = g.String(
                "Unknown D is the tangent at each point, kept in the right handle as an offset from the point while solving and starting at zero. C2 continuity gives h_next D[i-1] + 2 (h_prev + h_next) D[i] + h_prev D[i+1] = 3 (h_next / h_prev (K[i] - K[i-1]) + h_prev / h_next (K[i+1] - K[i])). Open ends have zero curvature: 2 D[0] + D[1] = 3 (K[1] - K[0]) / h and the mirror at the last point. The system is diagonally dominant, so each Jacobi step at least halves the error."
            )
            vector_math_2 = g.CurveHandlePositions().o.right - position
            evaluate_at_index_2 = vector_math_2.point.at(
                offset_point_in_curve.o.point_index
            )
            evaluate_at_index_3 = vector_math_2.point.at(
                offset_point_in_curve_1.o.point_index
            )
            vector_math_3 = (position.o.position - evaluate_at_index) * (
                switch_1 / switch
            ) + (evaluate_at_index_1 - position) * (switch / switch_1)
            vector_math_4 = (
                vector_math_3 * 3.0
                - evaluate_at_index_2 * switch_1
                - evaluate_at_index_3 * switch
            ) / ((switch + switch_1) * 2.0)
            vector_math_5 = (
                (evaluate_at_index_1 - position) * 3.0 / switch_1 - evaluate_at_index_3
            ) / (2.0, 2.0, 2.0)
            vector_math_6 = (
                (position.o.position - evaluate_at_index) * 3.0 / switch
                - evaluate_at_index_2
            ) / (2.0, 2.0, 2.0)
            switch_2 = boolean_math.switch.vector(
                boolean_math_1.switch.vector(vector_math_4, vector_math_6),
                vector_math_5,
            )
            (
                g.SetHandlePositions.right(
                    geometry.current, compare, position.o.position + switch_2
                )
                >> geometry.next
            )
        with g.Frame("Handles from tangents"):
            _string_2 = g.String(
                "A cubic Bézier segment of parameter length h leaves its start along h D / 3, so the right handle is K + h_next D / 3 and the left is K - h_prev D / 3."
            )
            vector_math_7 = g.CurveHandlePositions().o.right - position
            (
                geometry.result
                >> g.SetHandlePositions(
                    selection=compare,
                    position=position.o.position
                    - vector_math_7 * switch / (3.0, 3.0, 3.0),
                )
                >> g.SetHandlePositions.right(
                    selection=compare,
                    position=position.o.position
                    + vector_math_7 * switch_1 / (3.0, 3.0, 3.0),
                )
                >> curve_1
            )


ASSET = CurveSplineC2

ASSET_METADATA = {
    "description": "Convert curves to Bézier curves that pass through every control point with continuous curvature (a C2 cubic spline). Open curves have zero curvature at their ends. Curves with fewer than 3 points are converted to Bézier but not solved. Chord Length spaces the spline by the distance between points. Port of SolveCubicC2 from github.com/aurimasg/splines (MIT, Aurimas Gasiulis)",
    "catalog_id": "9c167a5c-d0a6-457d-9e9c-90f3edd29e10",
}
