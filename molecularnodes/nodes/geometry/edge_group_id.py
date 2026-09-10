# Node-group asset "Edge Group ID" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputInteger


class EdgeGroupID(AssetGeometryGroup):
    """
    Edge Group ID

    Parameters
    ----------
    group_id : InputInteger
        Group ID

    Inputs
    ------
    i.group_id : IntegerSocket
        Group ID

    Outputs
    -------
    o.difference : IntegerSocket
        Difference
    o.is_equal : BooleanSocket
        Is Equal
    """

    _name = "Edge Group ID"
    _asset_name = "Edge Group ID"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        group_id: IntegerSocket
        """Group ID"""

    class _Outputs(SocketAccessor):
        difference: IntegerSocket
        """Difference"""
        is_equal: BooleanSocket
        """Is Equal"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        group_id: InputInteger = 0,
    ):
        super().__init__(**{"Group ID": group_id})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        group_id = tree.inputs.integer("Group ID", 0, hide_value=True)
        difference = tree.outputs.integer("Difference", attribute_domain="EDGE")
        is_equal = tree.outputs.boolean("Is Equal", attribute_domain="EDGE")

        edge_vertices = g.EdgeVertices()
        evaluate_at_index = group_id.point.at(edge_vertices.o.vertex_index_1)
        evaluate_at_index_1 = group_id.point.at(edge_vertices.o.vertex_index_2)
        g.Compare.integer.equal(evaluate_at_index, evaluate_at_index_1) >> is_equal
        abs(evaluate_at_index - evaluate_at_index_1) >> difference


ASSET = EdgeGroupID

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
