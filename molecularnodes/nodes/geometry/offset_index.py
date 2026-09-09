# Node-group asset 'Offset Index' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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


class OffsetIndex(AssetGeometryGroup):
    """
    Offset Index

    Parameters
    ----------
    index : InputInteger
        The `Index` at which to evaluate this offset from
    offset : InputInteger
        The `Offset` to apply to the `Index` of the point

    Inputs
    ------
    i.index : IntegerSocket
        The `Index` at which to evaluate this offset from
    i.offset : IntegerSocket
        The `Offset` to apply to the `Index` of the point

    Outputs
    -------
    o.index : IntegerSocket
        The `Index` + `Offset`
    """

    _name = "Offset Index"
    _asset_name = "Offset Index"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.offset_index"}

    class _Inputs(SocketAccessor):
        index: IntegerSocket
        """The `Index` at which to evaluate this offset from"""
        offset: IntegerSocket
        """The `Offset` to apply to the `Index` of the point"""

    class _Outputs(SocketAccessor):
        index: IntegerSocket
        """The `Index` + `Offset`"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        index: InputInteger = 0,
        offset: InputInteger = 0,
    ):
        super().__init__(**{"Index": index, "Offset": offset})

    def _build_group(self, tree):
        index = tree.inputs.integer(
            "Index",
            0,
            description="The `Index` at which to evaluate this offset from",
            min_value=0,
            default_input="INDEX",
        )
        offset = tree.inputs.integer(
            "Offset", 0, description="The `Offset` to apply to the `Index` of the point"
        )
        index_1 = tree.outputs.integer("Index", description="The `Index` + `Offset`")

        index + offset >> index_1


ASSET = OffsetIndex

ASSET_METADATA = {
    "catalog_id": "7ccb8802-a69f-483e-bf6e-4a47aaa9e940",
}
