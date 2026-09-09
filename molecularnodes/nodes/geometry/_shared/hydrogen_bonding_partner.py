# Node group 'Hydrogen Bonding Partner' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import CustomGeometryGroup, IntegerSocket, SocketAccessor
from ..is_hydrogen import IsHydrogen


class HydrogenBondingPartner(CustomGeometryGroup):
    """
    Hydrogen Bonding Partner

    Outputs
    -------
    o.index : IntegerSocket
        Index
    """

    _name = "Hydrogen Bonding Partner"
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        index: IntegerSocket
        """Index"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self):
        super().__init__()

    def _build_group(self, tree):
        index = tree.outputs.integer("Index")

        edge_vertices = g.EdgeVertices()
        edges_of_vertex = g.EdgesOfVertex()
        evaluate_at_index = edge_vertices.o.vertex_index_1.edge.at(
            edges_of_vertex.o.edge_index
        )
        switch = g.Compare.integer.equal(
            evaluate_at_index, g.Index()
        ).o.result.switch.integer(
            evaluate_at_index,
            edge_vertices.o.vertex_index_2.edge.at(edges_of_vertex.o.edge_index),
        )
        (
            IsHydrogen(
                and_=g.Compare.integer.equal(g.EdgesOfVertex().o.total, 1)
            ).o.selection.switch.integer(g.Index(), switch)
            >> index
        )
