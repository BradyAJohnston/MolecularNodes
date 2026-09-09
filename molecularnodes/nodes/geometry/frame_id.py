# Node-group asset 'Frame ID' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy.builder import (
    AssetGeometryGroup,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputInteger
from ._shared.attribute_at_index import AttributeAtIndex


class FrameID(AssetGeometryGroup):
    """
    Read the `frame_id` attribute, created when multiple frames from a trajectory are merged into a single structure

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
    o.frame_id : IntegerSocket
        Read the `frame_id` attribute from the geometry
    """

    _name = "Frame ID"
    _asset_name = "Frame ID"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {
        "description": "Read the `frame_id` attribute, created when multiple frames from a trajectory are merged into a single structure",
        "node_tool_idname": "geometry.residue_id",
    }

    class _Inputs(SocketAccessor):
        index: IntegerSocket
        """Index"""

    class _Outputs(SocketAccessor):
        frame_id: IntegerSocket
        """Read the `frame_id` attribute from the geometry"""

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
        index = tree.inputs.integer(
            "Index", 0, min_value=0, hide_value=True, default_input="INDEX"
        )
        frame_id = tree.outputs.integer(
            "frame_id", description="Read the `frame_id` attribute from the geometry"
        )

        AttributeAtIndex(index=index, name="frame_id") >> frame_id


ASSET = FrameID

ASSET_METADATA = {
    "description": "Read the `frame_id` attribute, created when multiple frames from a trajectory are merged into a single structure",
    "catalog_id": "dfef0d3c-e718-420a-8b22-e7c3a3a9e333",
}
