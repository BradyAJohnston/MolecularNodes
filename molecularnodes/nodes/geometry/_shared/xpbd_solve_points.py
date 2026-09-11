# Node group "XPBD Solve Points" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
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
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputFloat, InputGeometry
from ..residue_id import ResidueID
from .constraint_distance import ConstraintDistance
from .inverse_mass import InverseMass


class XPBDSolvePoints(CustomGeometryGroup):
    """
    XPBD Solve Points

    Parameters
    ----------
    geometry : InputGeometry
        Geometry
    selection : InputBoolean
        Selection
    radius : InputFloat
        Radius
    alpha : InputFloat
        alpha
    deltat : InputFloat
        deltaT

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry
    i.selection : BooleanSocket
        Selection
    i.radius : FloatSocket
        Radius
    i.alpha : FloatSocket
        alpha
    i.deltat : FloatSocket
        deltaT

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "XPBD Solve Points"
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""
        selection: BooleanSocket
        """Selection"""
        radius: FloatSocket
        """Radius"""
        alpha: FloatSocket
        deltat: FloatSocket
        """deltaT"""

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
        geometry: InputGeometry = None,
        selection: InputBoolean = True,
        radius: InputFloat = 0.0,
        alpha: InputFloat = 0.0,
        deltat: InputFloat = 0.0,
    ):
        super().__init__(
            **{
                "Geometry": geometry,
                "Selection": selection,
                "Radius": radius,
                "alpha": alpha,
                "deltaT": deltat,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry("Geometry")
        selection = tree.inputs.boolean("Selection", True, hide_value=True)
        radius = tree.inputs.float("Radius", 0.0)
        alpha = tree.inputs.float("alpha", 0.0, min_value=-10_000.0, max_value=10_000.0)
        deltat = tree.inputs.float(
            "deltaT",
            0.0,
            min_value=-10_000.0,
            max_value=10_000.0,
            structure_type="SINGLE",
            force_non_field=True,
        )
        geometry_1 = tree.outputs.geometry("Geometry")

        group = ResidueID()
        group_1 = InverseMass()
        vector_math = (
            g.Position().o.position
            + g.RandomValue.vector((-1.0, -1.0, -1.0), (1.0, 1.0, 1.0)).o.value * 0.001
        )
        capture = g.CaptureAttribute.point(geometry=geometry)
        capture.items.vector("Vector", vector_math)
        capture_1 = g.CaptureAttribute.point(geometry=capture.o.geometry)
        index = capture_1.items.integer("Index", g.IndexOfNearest().o.index)
        evaluate_at_index = group.o.res_id.point.at(index.output)
        _compare = g.Compare.integer.not_equal(
            abs(group.o.res_id - evaluate_at_index), 1
        )
        math_1 = radius.point.at(index.output) + radius
        group_2 = ConstraintDistance(
            target=g.Position().o.position.point.at(index.output),
            distance=math_1,
            w1=group_1.o.w,
            w2=group_1.o.w.point.at(index.output),
            alpha=alpha,
            deltat=deltat,
        )
        compare_1 = g.Compare.integer.not_equal(
            g.ShortestEdgePaths(
                end_vertex=g.Compare.integer.equal(evaluate_at_index, g.Index())
            ).o.next_vertex_index,
            index.output,
        )
        (
            capture_1.o.geometry
            >> g.SetPosition(
                selection=compare_1.o.result & (selection & (math_1 > group_2.o.value)),
                offset=group_2.o.correction,
            )
            >> geometry_1
        )
