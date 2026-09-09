# Node-group asset 'Assembly ID' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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


class AssemblyID(AssetGeometryGroup):
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
    o.assembly_id : IntegerSocket
        Read the `frame_id` attribute from the geometry
    """

    _name = "Assembly ID"
    _asset_name = "Assembly ID"
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
        assembly_id: IntegerSocket
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
        assembly_id = tree.outputs.integer(
            "assembly_id", description="Read the `frame_id` attribute from the geometry"
        )

        AttributeAtIndex(index=index, name="assembly_id") >> assembly_id


ASSET = AssemblyID

ASSET_METADATA = {
    "description": "Read the `assembly_id` attribute, which defines which biological assembly the points belong to",
}
