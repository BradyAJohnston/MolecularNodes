# Node group "Curve Control Points" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
    GeometrySocket,
    IntegerSocket,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputGeometry


class CurveControlPoints(CustomGeometryGroup):
    """
    Curve Control Points

    Parameters
    ----------
    curve : InputGeometry
        Curve

    Inputs
    ------
    i.curve : GeometrySocket
        Curve

    Outputs
    -------
    o.points : GeometrySocket
        Control points in curve order, carrying their attributes
    o.handle_left : VectorSocket
        Handle Left
    o.handle_right : VectorSocket
        Handle Right
    o.curve_index : IntegerSocket
        Curve Index
    o.cyclic : BooleanSocket
        Cyclic
    """

    _name = "Curve Control Points"
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        curve: GeometrySocket
        """Curve"""

    class _Outputs(SocketAccessor):
        points: GeometrySocket
        """Control points in curve order, carrying their attributes"""
        handle_left: VectorSocket
        """Handle Left"""
        handle_right: VectorSocket
        """Handle Right"""
        curve_index: IntegerSocket
        """Curve Index"""
        cyclic: BooleanSocket
        """Cyclic"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        curve: InputGeometry = None,
    ):
        super().__init__(**{"Curve": curve})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        curve = tree.inputs.geometry("Curve")
        points = tree.outputs.geometry(
            "Points",
            description="Control points in curve order, carrying their attributes",
        )
        handle_left = tree.outputs.vector("Handle Left")
        handle_right = tree.outputs.vector("Handle Right")
        curve_index = tree.outputs.integer("Curve Index")
        cyclic = tree.outputs.boolean("Cyclic")

        capture = g.CaptureAttribute.point(geometry=g.SetSplineType.bezier(curve))
        handle_left_1 = capture.items.vector(
            g.CurveHandlePositions().o.left, "Handle Left"
        )
        handle_right_1 = capture.items.vector(
            g.CurveHandlePositions().o.right, "Handle Right"
        )
        curve_index_1 = capture.items.integer(
            g.CurveOfPoint().o.curve_index, "Curve Index"
        )
        cyclic_1 = capture.items.boolean(g.IsSplineCyclic(), "Cyclic")
        capture.o.geometry >> g.SetSplineType() >> g.CurveToPoints.evaluated() >> points

        handle_left_1.output >> handle_left
        handle_right_1.output >> handle_right
        curve_index_1.output >> curve_index
        cyclic_1.output >> cyclic
