# Node-group asset "Residue ID" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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


class ResidueID(AssetGeometryGroup):
    """
    Residue ID

    Outputs
    -------
    o.res_id : IntegerSocket
        Read the `res_id` attribute from the geometry
    """

    _name = "Residue ID"
    _asset_name = "Residue ID"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.residue_id"}

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        res_id: IntegerSocket
        """Read the `res_id` attribute from the geometry"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self) -> None:
        super().__init__()

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        res_id = tree.outputs.integer(
            "res_id", description="Read the `res_id` attribute from the geometry"
        )

        named_attribute = g.NamedAttribute.integer("res_id")

        named_attribute >> res_id


ASSET = ResidueID

ASSET_METADATA = {
    "catalog_id": "dfef0d3c-e718-420a-8b22-e7c3a3a9e333",
}
