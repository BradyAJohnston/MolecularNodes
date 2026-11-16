# Node-group asset "Curve Bridge" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
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


class CurveBridge(AssetGeometryGroup):
    """
    Join each selected curve to the next curve by index with a Bézier segment that leaves the end of one and arrives at the start of the next along their tangents. Point attributes are kept. Output curves are Bézier with free handles. Bridge Curves idea from github.com/Eterea/etr_curves_toolkit (GPL-3.0-or-later, Cristobal Vila)

    Parameters
    ----------
    curve : InputGeometry
        Curve
    selection : InputBoolean
        Curves to join to the next curve, evaluated on the curve domain. Cyclic curves are never joined
    tension : InputFloat
        Length of the bridge handles as a share of a third of the gap. 0 gives a straight bridge

    Inputs
    ------
    i.curve : GeometrySocket
        Curve
    i.selection : BooleanSocket
        Curves to join to the next curve, evaluated on the curve domain. Cyclic curves are never joined
    i.tension : FloatSocket
        Length of the bridge handles as a share of a third of the gap. 0 gives a straight bridge

    Outputs
    -------
    o.curve : GeometrySocket
        Curve
    o.bridge : BooleanSocket
        Points at either end of a new bridge segment
    """

    _name = "Curve Bridge"
    _asset_name = "Curve Bridge"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "description": "Join each selected curve to the next curve by index with a Bézier segment that leaves the end of one and arrives at the start of the next along their tangents. Point attributes are kept. Output curves are Bézier with free handles. Bridge Curves idea from github.com/Eterea/etr_curves_toolkit (GPL-3.0-or-later, Cristobal Vila)"
    }

    class _Inputs(SocketAccessor):
        curve: GeometrySocket
        """Curve"""
        selection: BooleanSocket
        """Curves to join to the next curve, evaluated on the curve domain. Cyclic curves are never joined"""
        tension: FloatSocket
        """Length of the bridge handles as a share of a third of the gap. 0 gives a straight bridge"""

    class _Outputs(SocketAccessor):
        curve: GeometrySocket
        """Curve"""
        bridge: BooleanSocket
        """Points at either end of a new bridge segment"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        curve: InputGeometry = None,
        selection: InputBoolean = True,
        tension: InputFloat = 1.0,
    ):
        super().__init__(**{"Curve": curve, "Selection": selection, "Tension": tension})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        curve = tree.inputs.geometry("Curve")
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Curves to join to the next curve, evaluated on the curve domain. Cyclic curves are never joined",
            hide_value=True,
        )
        tension = tree.inputs.float(
            "Tension",
            1.0,
            description="Length of the bridge handles as a share of a third of the gap. 0 gives a straight bridge",
            min_value=0.0,
        )
        curve_1 = tree.outputs.geometry("Curve")
        bridge = tree.outputs.boolean(
            "Bridge", description="Points at either end of a new bridge segment"
        )

        set_spline_type = g.SetSplineType.bezier(curve)
        with g.Frame("Which curves join"):
            index = g.Index()
            is_spline_cyclic = g.IsSplineCyclic()
            boolean_math = (
                selection
                & ~is_spline_cyclic.o.cyclic
                & ~is_spline_cyclic.o.cyclic.spline.at(index.o.index + 1)
            )
            boolean_math_1 = boolean_math & (
                index
                < g.DomainSize(
                    geometry=set_spline_type, component="CURVE"
                ).o.spline_count
                - 1
            )
            _string = g.String(
                "A curve joins the next one when it is selected and neither is cyclic. Each run of joined curves shares one group ID: the count of unjoined curves before it."
            )
            capture = g.CaptureAttribute.curve(geometry=set_spline_type)
            join = capture.items.boolean(boolean_math_1, "Join")
            group = capture.items.integer(
                boolean_math_1.switch.integer(1).spline.trailing(), "Group"
            )
        with g.Frame("Bridge handles"):
            curve_tangent = g.CurveTangent()
            curve_of_point = g.CurveOfPoint()
            compare = g.Compare.integer.equal(
                curve_of_point.o.index_in_curve,
                g.PointsOfCurve(curve_of_point.o.curve_index).o.total - 1,
            )
            position = g.Position()
            _string_1 = g.String(
                "The end of a joined curve and the start of the next get free handles along their tangents, Tension times a third of the gap long, so the segment between them is the bridge."
            )
            boolean_math_2 = join.output & compare
            boolean_math_3 = join.output.spline.at(
                curve_of_point.o.curve_index - 1
            ) & g.Compare.integer.equal(curve_of_point.o.index_in_curve, 0)
            boolean_math_4 = boolean_math_3 & (curve_of_point.o.curve_index > 0)
            math_1 = (
                (
                    position.o.position.point.at(g.Index().o.index + 1) - position
                ).length()
                * tension
                / 3.0
            )
            math_2 = (
                position.o.position
                - position.o.position.point.at(g.Index().o.index - 1)
            ).length() * tension
            set_handle_positions = (
                g.SetHandleType.free(
                    capture.o.geometry, boolean_math_2 | boolean_math_4
                )
                >> g.SetHandlePositions.right(
                    selection=boolean_math_2,
                    position=position.o.position + curve_tangent.o.tangent * math_1,
                )
                >> g.SetHandlePositions(
                    selection=boolean_math_4,
                    position=position.o.position
                    - curve_tangent.o.tangent * (math_2 / 3.0),
                )
            )
            capture_1 = g.CaptureAttribute.point(geometry=set_handle_positions)
            bridge_1 = capture_1.items.boolean(
                boolean_math_2 | boolean_math_4, "Bridge"
            )
            group_1 = capture_1.items.integer(group.output, "Group")
        curve_control_points = CurveControlPoints(curve=capture_1.o.geometry)
        (
            ControlPointsToCurve(
                points=curve_control_points,
                curve_group_id=group_1.output,
                weight=g.Index(),
                handle_left=curve_control_points.o.handle_left,
                handle_right=curve_control_points.o.handle_right,
                cyclic=curve_control_points.o.cyclic,
            )
            >> curve_1
        )

        bridge_1.output >> bridge


ASSET = CurveBridge

ASSET_METADATA = {
    "description": "Join each selected curve to the next curve by index with a Bézier segment that leaves the end of one and arrives at the start of the next along their tangents. Point attributes are kept. Output curves are Bézier with free handles. Bridge Curves idea from github.com/Eterea/etr_curves_toolkit (GPL-3.0-or-later, Cristobal Vila)",
    "catalog_id": "9c167a5c-d0a6-457d-9e9c-90f3edd29e10",
}
