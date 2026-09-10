# Node-group asset "Edge Info" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputInteger


class EdgeInfo(AssetGeometryGroup):
    """
    Edge Info

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
    o.is_valid : BooleanSocket
        Whether there is a valid edge corresponding to the given index
    o.point_index : IntegerSocket
        The index for the other point involved in this edge, -1 if not connected
    o.point_position : VectorSocket
        The position for the other point involved in this edge, (0, 0, 0) if not connected
    o.edge_index : IntegerSocket
        The index on the edge domain for the selected edge. -1 if not connected
    o.edge_vector : VectorSocket
        The vector along the selected edge. (0, 0, 0) if not connected
    o.edge_length : FloatSocket
        Length of the selected edge, -1 if not connected
    """

    _name = "Edge Info"
    _asset_name = "Edge Info"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.edge_info"}

    class _Inputs(SocketAccessor):
        vertex_index: IntegerSocket
        """Vertex Index"""
        edge_index: IntegerSocket
        """Index within the gorup of edges that are connected to this point"""

    class _Outputs(SocketAccessor):
        is_valid: BooleanSocket
        """Whether there is a valid edge corresponding to the given index"""
        point_index: IntegerSocket
        """The index for the other point involved in this edge, -1 if not connected"""
        point_position: VectorSocket
        """The position for the other point involved in this edge, (0, 0, 0) if not connected"""
        edge_index: IntegerSocket
        """The index on the edge domain for the selected edge. -1 if not connected"""
        edge_vector: VectorSocket
        """The vector along the selected edge. (0, 0, 0) if not connected"""
        edge_length: FloatSocket
        """Length of the selected edge, -1 if not connected"""

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
        is_valid = tree.outputs.boolean(
            "Is Valid",
            description="Whether there is a valid edge corresponding to the given index",
        )
        point_index = tree.outputs.integer(
            "Point Index",
            -1,
            description="The index for the other point involved in this edge, -1 if not connected",
            min_value=-1,
        )
        point_position = tree.outputs.vector(
            "Point Position",
            description="The position for the other point involved in this edge, (0, 0, 0) if not connected",
        )
        edge_index_1 = tree.outputs.integer(
            "Edge Index",
            -1,
            description="The index on the edge domain for the selected edge. -1 if not connected",
            min_value=-1,
        )
        edge_vector = tree.outputs.vector(
            "Edge Vector",
            description="The vector along the selected edge. (0, 0, 0) if not connected",
            subtype="EULER",
        )
        edge_length = tree.outputs.float(
            "Edge Length",
            -1.0,
            description="Length of the selected edge, -1 if not connected",
            min_value=0.0,
        )

        edge_vertices = g.EdgeVertices()
        edges_of_vertex = g.EdgesOfVertex(
            vertex_index=vertex_index, sort_index=edge_index
        )
        with g.Frame("Check edge exists, or return index -1 and (0, 0, 0) vector"):
            compare = edge_index < edges_of_vertex.o.total
        evaluate_at_index = g.EvaluateAtIndex.edge.integer(
            g.Math(
                value=edge_vertices.o.vertex_index_1,
                value_001=edge_vertices.o.vertex_index_2,
            ),
            edges_of_vertex.o.edge_index,
        )
        math_1 = g.Math.subtract(evaluate_at_index, g.Index())
        position = g.Position()
        evaluate_at_index_1 = position.o.position.point.at(math_1)
        compare.switch.integer(-1, math_1) >> point_index
        compare.switch.integer(-1, edges_of_vertex.o.edge_index) >> edge_index_1
        vector_math = evaluate_at_index_1 - position
        compare.switch.vector((0.0, 0.0, 0.0), vector_math) >> edge_vector
        compare.switch.float(-1.0, vector_math.length()) >> edge_length

        compare >> is_valid
        evaluate_at_index_1 >> point_position


ASSET = EdgeInfo

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
