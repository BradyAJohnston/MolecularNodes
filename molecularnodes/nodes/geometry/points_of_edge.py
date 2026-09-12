# Node-group asset "Points of Edge" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputInteger
from .edge_info import EdgeInfo


class PointsOfEdge(AssetGeometryGroup):
    """
    Points of Edge

    Parameters
    ----------
    vertex_index : InputInteger
        Vertex Index
    edge_index : InputInteger
        Index within the gorup of edges that are connected to this point

    Inputs
    ------
    i.vertex_index : IntegerSocket
        Vertex Index
    i.edge_index : IntegerSocket
        Index within the gorup of edges that are connected to this point

    Outputs
    -------
    o._0 : IntegerSocket
        Index for the 0th point, connected to the point at the end of the selected edge. Returns -1 if not connected or self
    o._1 : IntegerSocket
        Index for the 1th point, connected to the point at the end of the selected edge. Returns -1 if not connected or self
    o._2 : IntegerSocket
        Index for the 2th point, connected to the point at the end of the selected edge. Returns -1 if not connected or self
    o._3 : IntegerSocket
        Index for the 3th point, connected to the point at the end of the selected edge. Returns -1 if not connected or self
    o.total : IntegerSocket
        Number of edges conncted to the connected point, including this edge
    """

    _name = "Points of Edge"
    _asset_name = "Points of Edge"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.points_of_edge"}

    class _Inputs(SocketAccessor):
        vertex_index: IntegerSocket
        """Vertex Index"""
        edge_index: IntegerSocket
        """Index within the gorup of edges that are connected to this point"""

    class _Outputs(SocketAccessor):
        _0: IntegerSocket
        """Index for the 0th point, connected to the point at the end of the selected edge. Returns -1 if not connected or self"""
        _1: IntegerSocket
        """Index for the 1th point, connected to the point at the end of the selected edge. Returns -1 if not connected or self"""
        _2: IntegerSocket
        """Index for the 2th point, connected to the point at the end of the selected edge. Returns -1 if not connected or self"""
        _3: IntegerSocket
        """Index for the 3th point, connected to the point at the end of the selected edge. Returns -1 if not connected or self"""
        total: IntegerSocket
        """Number of edges conncted to the connected point, including this edge"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        vertex_index: InputInteger = 0,
        edge_index: InputInteger = 0,
    ):
        super().__init__(**{"Vertex Index": vertex_index, "Edge Index": edge_index})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        vertex_index = tree.inputs.integer(
            "Vertex Index", 0, hide_value=True, default_input="INDEX"
        )
        edge_index = tree.inputs.integer(
            "Edge Index",
            0,
            description="Index within the gorup of edges that are connected to this point",
            min_value=0,
            max_value=3,
        )
        n_0 = tree.outputs.integer(
            "0",
            -1,
            description="Index for the 0th point, connected to the point at the end of the selected edge. Returns -1 if not connected or self",
            min_value=-1,
        )
        n_1 = tree.outputs.integer(
            "1",
            -1,
            description="Index for the 1th point, connected to the point at the end of the selected edge. Returns -1 if not connected or self",
            min_value=-1,
        )
        n_2 = tree.outputs.integer(
            "2",
            -1,
            description="Index for the 2th point, connected to the point at the end of the selected edge. Returns -1 if not connected or self",
            min_value=-1,
        )
        n_3 = tree.outputs.integer(
            "3",
            -1,
            description="Index for the 3th point, connected to the point at the end of the selected edge. Returns -1 if not connected or self",
            min_value=-1,
        )
        total = tree.outputs.integer(
            "Total",
            description="Number of edges conncted to the connected point, including this edge",
            min_value=0,
        )

        group = EdgeInfo(vertex_index=vertex_index, edge_index=edge_index)
        g.EdgesOfVertex().o.total.point.at(group.o.point_index) >> total
        evaluate_at_index = EdgeInfo(vertex_index=vertex_index).o.point_index.point.at(
            group.o.point_index
        )
        evaluate_at_index_1 = EdgeInfo(
            vertex_index=vertex_index, edge_index=1
        ).o.point_index.point.at(group.o.point_index)
        evaluate_at_index_2 = EdgeInfo(
            vertex_index=vertex_index, edge_index=2
        ).o.point_index.point.at(group.o.point_index)
        evaluate_at_index_3 = EdgeInfo(
            vertex_index=edge_index, edge_index=3
        ).o.point_index.point.at(group.o.point_index)
        index = g.Index()
        with g.Frame("check if selecting self, return -1 if so"):
            (
                g.Compare.integer.equal(
                    evaluate_at_index, index
                ).o.result.switch.integer(evaluate_at_index, -1)
                >> n_0
            )
            (
                g.Compare.integer.equal(
                    evaluate_at_index_1, index
                ).o.result.switch.integer(evaluate_at_index_1, -1)
                >> n_1
            )
            (
                g.Compare.integer.equal(
                    evaluate_at_index_2, index
                ).o.result.switch.integer(evaluate_at_index_2, -1)
                >> n_2
            )
            (
                g.Compare.integer.equal(
                    evaluate_at_index_3, index
                ).o.result.switch.integer(evaluate_at_index_3, -1)
                >> n_3
            )


ASSET = PointsOfEdge

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
