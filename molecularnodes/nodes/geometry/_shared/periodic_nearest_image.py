# Node group "Periodic Nearest Image" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import CustomGeometryGroup, SocketAccessor, VectorSocket
from nodebpy.types import InputVector


class PeriodicNearestImage(CustomGeometryGroup):
    """
    Periodic Nearest Image

    Parameters
    ----------
    vector : InputVector
        Vectors to find the shortest periodic image of
    a : InputVector
        Box vector A
    b : InputVector
        Box vector B
    c : InputVector
        Box vector C

    Inputs
    ------
    i.vector : VectorSocket
        Vectors to find the shortest periodic image of
    i.a : VectorSocket
        Box vector A
    i.b : VectorSocket
        Box vector B
    i.c : VectorSocket
        Box vector C

    Outputs
    -------
    o.vector : VectorSocket
        The shortest periodic image of each vector
    """

    _name = "Periodic Nearest Image"
    _color_tag = "VECTOR"

    class _Inputs(SocketAccessor):
        vector: VectorSocket
        """Vectors to find the shortest periodic image of"""
        a: VectorSocket
        """Box vector A"""
        b: VectorSocket
        """Box vector B"""
        c: VectorSocket
        """Box vector C"""

    class _Outputs(SocketAccessor):
        vector: VectorSocket
        """The shortest periodic image of each vector"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        vector: InputVector = None,
        a: InputVector = None,
        b: InputVector = None,
        c: InputVector = None,
    ):
        super().__init__(**{"Vector": vector, "A": a, "B": b, "C": c})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        vector = tree.inputs.vector(
            "Vector",
            (0.0, 0.0, 0.0),
            description="Vectors to find the shortest periodic image of",
            hide_value=True,
        )
        a = tree.inputs.vector("A", (0.0, 0.0, 0.0), description="Box vector A")
        b = tree.inputs.vector("B", (0.0, 0.0, 0.0), description="Box vector B")
        c_ = tree.inputs.vector("C", (0.0, 0.0, 0.0), description="Box vector C")
        vector_1 = tree.outputs.vector(
            "Vector", description="The shortest periodic image of each vector"
        )

        with g.Frame("Lattice Points"):
            index = g.Index()
            _string = g.String(
                "The 27 lattice points of shifts -1, 0 and 1 along A, B and C, numbered by index."
            )
            vector_math = (
                a * (index.o.index % 3 - 1)
                + b * (index.o.index // 3 % 3 - 1)
                + c_ * (index.o.index // 9 - 1)
            )
            points = g.Points(position=vector_math, count=27)
        with g.Frame("Staged Rounding"):
            _string_1 = g.String(
                "Remove whole box vectors C, then B, then A, as the box matrix is lower triangular. This leaves each vector within about one cell of its shortest image."
            )
            vector_math_1 = vector - c_ * (vector.z / c_.z).round()
            vector_math_2 = vector_math_1 - b * (vector_math_1.y / b.y).round()
            vector_math_3 = vector_math_2 - a * (vector_math_2.x / a.x).round()
        with g.Frame("Nearest Lattice Point"):
            _string_2 = g.String(
                "The nearest lattice point to the vector, found with a KD-tree, rebuilt from its index, is the shift that gives the shortest image."
            )
            sample_nearest = g.SampleNearest.point(points, vector_math_3)
            vector_math_4 = (
                a * (sample_nearest.o.index % 3 - 1)
                + b * (sample_nearest.o.index // 3 % 3 - 1)
                + c_ * (sample_nearest.o.index // 9 - 1)
            )
        vector_math_3 - vector_math_4 >> vector_1
