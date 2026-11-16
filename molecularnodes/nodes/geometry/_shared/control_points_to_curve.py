# Node group "Control Points to Curve" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    BooleanSocket,
    CustomGeometryGroup,
    FloatSocket,
    GeometrySocket,
    IntegerSocket,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import (
    InputBoolean,
    InputFloat,
    InputGeometry,
    InputInteger,
    InputVector,
)


class ControlPointsToCurve(CustomGeometryGroup):
    """
    Control Points to Curve

    Parameters
    ----------
    points : InputGeometry
        Points
    curve_group_id : InputInteger
        Curve each point joins
    weight : InputFloat
        Order of points in a curve
    handle_left : InputVector
        Handle Left
    handle_right : InputVector
        Handle Right
    cyclic : InputBoolean
        Cyclic

    Inputs
    ------
    i.points : GeometrySocket
        Points
    i.curve_group_id : IntegerSocket
        Curve each point joins
    i.weight : FloatSocket
        Order of points in a curve
    i.handle_left : VectorSocket
        Handle Left
    i.handle_right : VectorSocket
        Handle Right
    i.cyclic : BooleanSocket
        Cyclic

    Outputs
    -------
    o.curve : GeometrySocket
        Curve
    """

    _name = "Control Points to Curve"
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        points: GeometrySocket
        """Points"""
        curve_group_id: IntegerSocket
        """Curve each point joins"""
        weight: FloatSocket
        """Order of points in a curve"""
        handle_left: VectorSocket
        """Handle Left"""
        handle_right: VectorSocket
        """Handle Right"""
        cyclic: BooleanSocket
        """Cyclic"""

    class _Outputs(SocketAccessor):
        curve: GeometrySocket
        """Curve"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        points: InputGeometry = None,
        curve_group_id: InputInteger = 0,
        weight: InputFloat = 0.0,
        handle_left: InputVector = None,
        handle_right: InputVector = None,
        cyclic: InputBoolean = False,
    ):
        super().__init__(
            **{
                "Points": points,
                "Curve Group ID": curve_group_id,
                "Weight": weight,
                "Handle Left": handle_left,
                "Handle Right": handle_right,
                "Cyclic": cyclic,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        points = tree.inputs.geometry("Points")
        curve_group_id = tree.inputs.integer(
            "Curve Group ID", 0, description="Curve each point joins", hide_value=True
        )
        weight = tree.inputs.float(
            "Weight", 0.0, description="Order of points in a curve", hide_value=True
        )
        handle_left = tree.inputs.vector(
            "Handle Left", (0.0, 0.0, 0.0), hide_value=True
        )
        handle_right = tree.inputs.vector(
            "Handle Right", (0.0, 0.0, 0.0), hide_value=True
        )
        cyclic = tree.inputs.boolean("Cyclic", False, hide_value=True)
        curve = tree.outputs.geometry("Curve")

        set_spline_type = (
            points
            >> g.PointsToCurves(curve_group_id=curve_group_id, weight=weight)
            >> g.SetSplineType.bezier()
        )
        (
            g.SetHandleType.free(set_spline_type)
            >> g.SetHandlePositions(position=handle_left)
            >> g.SetHandlePositions.right(position=handle_right)
            >> g.SetSplineCyclic(cyclic=cyclic)
            >> curve
        )
