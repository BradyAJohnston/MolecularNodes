# Node-group asset "Fragment ID" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
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


class FragmentID(AssetGeometryGroup):
    """
    Integer ID of the bonded molecule (MDAnalysis fragment) each atom belongs to, for structures with bonds

    Outputs
    -------
    o.fragment_id : IntegerSocket
        The `fragment_id` attribute read from the points, numbering the bonded molecules from 0
    """

    _name = "Fragment ID"
    _asset_name = "Fragment ID"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"
    _tree_properties = {
        "description": "Integer ID of the bonded molecule (MDAnalysis fragment) each atom belongs to, for structures with bonds",
        "node_tool_idname": "geometry.fragment_id",
    }

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        fragment_id: IntegerSocket
        """The `fragment_id` attribute read from the points, numbering the bonded molecules from 0"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self) -> None:
        super().__init__()

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        fragment_id = tree.outputs.integer(
            "fragment_id",
            description="The `fragment_id` attribute read from the points, numbering the bonded molecules from 0",
        )

        named_attribute = g.NamedAttribute.integer("fragment_id")

        named_attribute >> fragment_id


ASSET = FragmentID

ASSET_METADATA = {
    "description": "Integer ID of the bonded molecule (MDAnalysis fragment) each atom belongs to, for structures with bonds",
    "catalog_id": "dfef0d3c-e718-420a-8b22-e7c3a3a9e333",
}
