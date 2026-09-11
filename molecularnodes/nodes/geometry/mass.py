# Node-group asset "Mass" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputInteger
from .fallback_float import FallbackFloat


class Mass(AssetGeometryGroup):
    """
    Mass

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
    o.mass : FloatSocket
        Read the `mass` attribute from the geometry
    """

    _name = "Mass"
    _asset_name = "Mass"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.mass"}

    class _Inputs(SocketAccessor):
        index: IntegerSocket
        """Index"""

    class _Outputs(SocketAccessor):
        mass: FloatSocket
        """Read the `mass` attribute from the geometry"""

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

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        index = tree.inputs.integer("Index", 0, min_value=0, default_input="INDEX")
        mass = tree.outputs.float(
            "mass", description="Read the `mass` attribute from the geometry"
        )

        FallbackFloat(name="mass", fallback=1.0).o.value.point.at(index) >> mass


ASSET = Mass

ASSET_METADATA = {
    "catalog_id": "dfef0d3c-e718-420a-8b22-e7c3a3a9e333",
}
