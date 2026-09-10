# Node-group asset 'Bond Count' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputInteger


class BondCount(AssetGeometryGroup):
    """
    Bond Count

    Parameters
    ----------
    index : InputInteger
        Index

    Inputs
    ------
    i.index : IntegerSocket
        Index

    Outputs
    -------
    o.is_bonded : BooleanSocket
        If the point has an edge / bond to another point
    o.bonds : IntegerSocket
        The number of bonds or edges that a point has
    """

    _name = "Bond Count"
    _asset_name = "Bond Count"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.bond_count"}

    class _Inputs(SocketAccessor):
        index: IntegerSocket
        """Index"""

    class _Outputs(SocketAccessor):
        is_bonded: BooleanSocket
        """If the point has an edge / bond to another point"""
        bonds: IntegerSocket
        """The number of bonds or edges that a point has"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        index: InputInteger = 0,
    ):
        super().__init__(**{"Index": index})

    def _build_group(self, tree):
        index = tree.inputs.integer("Index", 0, min_value=0, default_input="INDEX")
        is_bonded = tree.outputs.boolean(
            "Is Bonded", description="If the point has an edge / bond to another point"
        )
        bonds = tree.outputs.integer(
            "Bonds", description="The number of bonds or edges that a point has"
        )

        edges_of_vertex = g.EdgesOfVertex(vertex_index=index)
        (edges_of_vertex.o.total > 0) >> is_bonded

        edges_of_vertex.o.total >> bonds


ASSET = BondCount

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
