# Node-group asset 'Is Boundary Edge' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean


class IsBoundaryEdge(AssetGeometryGroup):
    """
    Is Boundary Edge

    Parameters
    ----------
    mask : InputBoolean
        Mask

    Inputs
    ------
    i.mask : BooleanSocket
        Mask

    Outputs
    -------
    o.selection : BooleanSocket
        The calculated selection
    """

    _name = "Is Boundary Edge"
    _asset_name = "Is Boundary Edge"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.is_boundary_edge"}

    class _Inputs(SocketAccessor):
        mask: BooleanSocket
        """Mask"""

    class _Outputs(SocketAccessor):
        selection: BooleanSocket
        """The calculated selection"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        mask: InputBoolean = True,
    ):
        super().__init__(**{"Mask": mask})

    def _build_group(self, tree):
        mask = tree.inputs.boolean("Mask", True, hide_value=True)
        selection = tree.outputs.boolean(
            "Selection", description="The calculated selection"
        )

        edge_vertices = g.EdgeVertices()
        boolean_math = mask.point.at(edge_vertices.o.vertex_index_1) & mask.point.at(
            edge_vertices.o.vertex_index_2
        )
        (
            (
                g.Compare.integer.equal(g.EdgeNeighbors(), 1).o.result
                & boolean_math.edge.evaluate()
            )
            >> selection
        )


ASSET = IsBoundaryEdge

ASSET_METADATA = {
    "catalog_id": "bd1f205b-fea5-4700-b2c2-754f3321e969",
}
