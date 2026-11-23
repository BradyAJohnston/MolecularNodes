# Node-group asset "Residue Name" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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


class ResidueName(AssetGeometryGroup):
    """
    Residue Name

    Outputs
    -------
    o.res_name : IntegerSocket
        Read the `res_name` attribute from the geometry
    """

    _name = "Residue Name"
    _asset_name = "Residue Name"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.residue_name"}

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        res_name: IntegerSocket
        """Read the `res_name` attribute from the geometry"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self) -> None:
        super().__init__()

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        res_name = tree.outputs.integer(
            "res_name", description="Read the `res_name` attribute from the geometry"
        )

        named_attribute = g.NamedAttribute.integer("res_name")

        named_attribute >> res_name


ASSET = ResidueName

ASSET_METADATA = {
    "catalog_id": "dfef0d3c-e718-420a-8b22-e7c3a3a9e333",
}
