# Node-group asset 'Residue ID' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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


class ResidueID(AssetGeometryGroup):
    """
    Residue ID

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
    o.res_id : IntegerSocket
        Read the `res_id` attribute from the geometry
    """

    _name = "Residue ID"
    _asset_name = "Residue ID"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.residue_id"}

    class _Inputs(SocketAccessor):
        index: IntegerSocket
        """Index"""

    class _Outputs(SocketAccessor):
        res_id: IntegerSocket
        """Read the `res_id` attribute from the geometry"""

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
        res_id = tree.outputs.integer(
            "res_id", description="Read the `res_id` attribute from the geometry"
        )

        AttributeAtIndex(index=index) >> res_id


ASSET = ResidueID

ASSET_METADATA = {
    "catalog_id": "dfef0d3c-e718-420a-8b22-e7c3a3a9e333",
}
