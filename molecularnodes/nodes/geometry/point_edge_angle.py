# Node-group asset "Point Edge Angle" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from .edge_info import EdgeInfo
from .vector_angle import VectorAngle


class PointEdgeAngle(AssetGeometryGroup):
    """
    Point Edge Angle

    Parameters
    ----------
    vertex_index : InputInteger
        The index of the point at which to evaluate this node
    edge_a : InputInteger
        The index of the edges of this point to select
    edge_b : InputInteger
        The index of the edges of this point to select

    Inputs
    ------
    i.vertex_index : IntegerSocket
        The index of the point at which to evaluate this node
    i.edge_a : IntegerSocket
        The index of the edges of this point to select
    i.edge_b : IntegerSocket
        The index of the edges of this point to select

    Outputs
    -------
    o.is_valid : BooleanSocket
        Both of the selected edges are valid
    o.angle : FloatSocket
        Angle between the two selected edges in radians. Returns 0 if not valid
    o.edge_index_a : IntegerSocket
        Index for "Edge A" in the Edge domain of the geometry. Returns -1 if not valid
    o.edge_index_b : IntegerSocket
        Index for "Edge B" in the Edge domain of the geometry. Returns -1 if not valid
    o.edge_vector_a : VectorSocket
        Vector from the current point to the other point in Edge A. Returns (0, 0, 0) if not valid.
    o.edge_vector_b : VectorSocket
        Vector from the current point to the other point in Edge B. Returns (0, 0, 0) if not valid.
    """

    _name = "Point Edge Angle"
    _asset_name = "Point Edge Angle"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.point_edge_angle"}

    class _Inputs(SocketAccessor):
        vertex_index: IntegerSocket
        """The index of the point at which to evaluate this node"""
        edge_a: IntegerSocket
        """The index of the edges of this point to select"""
        edge_b: IntegerSocket
        """The index of the edges of this point to select"""

    class _Outputs(SocketAccessor):
        is_valid: BooleanSocket
        """Both of the selected edges are valid"""
        angle: FloatSocket
        """Angle between the two selected edges in radians. Returns 0 if not valid"""
        edge_index_a: IntegerSocket
        """Index for "Edge A" in the Edge domain of the geometry. Returns -1 if not valid"""
        edge_index_b: IntegerSocket
        """Index for "Edge B" in the Edge domain of the geometry. Returns -1 if not valid"""
        edge_vector_a: VectorSocket
        """Vector from the current point to the other point in Edge A. Returns (0, 0, 0) if not valid."""
        edge_vector_b: VectorSocket
        """Vector from the current point to the other point in Edge B. Returns (0, 0, 0) if not valid."""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        vertex_index: InputInteger = 0,
        edge_a: InputInteger = 0,
        edge_b: InputInteger = 1,
    ):
        super().__init__(
            **{"Vertex Index": vertex_index, "Edge A": edge_a, "Edge B": edge_b}
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        vertex_index = tree.inputs.integer(
            "Vertex Index",
            0,
            description="The index of the point at which to evaluate this node",
            hide_value=True,
            default_input="INDEX",
        )
        edge_a = tree.inputs.integer(
            "Edge A",
            0,
            description="The index of the edges of this point to select",
            min_value=0,
            max_value=3,
        )
        edge_b = tree.inputs.integer(
            "Edge B",
            1,
            description="The index of the edges of this point to select",
            min_value=0,
            max_value=3,
        )
        is_valid = tree.outputs.boolean(
            "Is Valid", description="Both of the selected edges are valid"
        )
        angle = tree.outputs.float(
            "Angle",
            -1.0,
            description="Angle between the two selected edges in radians. Returns 0 if not valid",
        )
        edge_index_a = tree.outputs.integer(
            "Edge Index A",
            -1,
            description='Index for "Edge A" in the Edge domain of the geometry. Returns -1 if not valid',
            min_value=0,
        )
        edge_index_b = tree.outputs.integer(
            "Edge Index B",
            -1,
            description='Index for "Edge B" in the Edge domain of the geometry. Returns -1 if not valid',
            min_value=0,
        )
        edge_vector_a = tree.outputs.vector(
            "Edge Vector A",
            description="Vector from the current point to the other point in Edge A. Returns (0, 0, 0) if not valid.",
        )
        edge_vector_b = tree.outputs.vector(
            "Edge Vector B",
            description="Vector from the current point to the other point in Edge B. Returns (0, 0, 0) if not valid.",
        )

        group = EdgeInfo(vertex_index=vertex_index, edge_index=edge_a)
        group_1 = EdgeInfo(vertex_index=vertex_index, edge_index=edge_b)
        boolean_math = g.Compare.integer.not_equal(edge_a, edge_b).o.result & (
            group.o.is_valid & group_1.o.is_valid
        )
        (
            boolean_math.switch.float(
                true=VectorAngle(a=group.o.edge_vector, b=group_1.o.edge_vector).o.angle
            )
            >> angle
        )

        boolean_math >> is_valid
        group.o.edge_index >> edge_index_a
        group_1.o.edge_index >> edge_index_b
        group.o.edge_vector >> edge_vector_a
        group_1.o.edge_vector >> edge_vector_b


ASSET = PointEdgeAngle

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
