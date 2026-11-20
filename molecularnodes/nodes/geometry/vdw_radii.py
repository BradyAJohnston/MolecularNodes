# Node-group asset "VDW Radii" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    PackageLibrary,
    SocketAccessor,
)
from ._shared.mn_world_scale import MN_world_scale
from .fallback_float import FallbackFloat


class VDWRadii(AssetGeometryGroup):
    """
    VDW Radii

    Outputs
    -------
    o.vdw_radii : FloatSocket
        Read the `vdw_radii` attribute from the geometry
    """

    _name = "VDW Radii"
    _asset_name = "VDW Radii"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.vdw_radii"}

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        vdw_radii: FloatSocket
        """Read the `vdw_radii` attribute from the geometry"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self) -> None:
        super().__init__()

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        vdw_radii = tree.outputs.float(
            "vdw_radii", description="Read the `vdw_radii` attribute from the geometry"
        )

        FallbackFloat(name="vdw_radii", fallback=MN_world_scale()) >> vdw_radii


ASSET = VDWRadii

ASSET_METADATA = {
    "catalog_id": "dfef0d3c-e718-420a-8b22-e7c3a3a9e333",
}
