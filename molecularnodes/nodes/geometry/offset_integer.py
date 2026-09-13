# Node-group asset "Offset Integer" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    AssetGeometryGroup,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputInteger
from .offset_index import OffsetIndex


class OffsetInteger(AssetGeometryGroup):
    """
    Offset Integer

    Parameters
    ----------
    integer : InputInteger
        The field to evaluate at the given `Index` + `Offset` on the point domain
    index : InputInteger
        The `Index` at which to evaluate this offset from
    offset : InputInteger
        The offset to apply to the `Index` before evaluating the input field

    Inputs
    ------
    i.integer : IntegerSocket
        The field to evaluate at the given `Index` + `Offset` on the point domain
    i.index : IntegerSocket
        The `Index` at which to evaluate this offset from
    i.offset : IntegerSocket
        The offset to apply to the `Index` before evaluating the input field

    Outputs
    -------
    o.integer : IntegerSocket
        The field evaluated at the offset `Index` value
    """

    _name = "Offset Integer"
    _asset_name = "Offset Integer"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {"node_tool_idname": "geometry.offset_integer"}

    class _Inputs(SocketAccessor):
        integer: IntegerSocket
        """The field to evaluate at the given `Index` + `Offset` on the point domain"""
        index: IntegerSocket
        """The `Index` at which to evaluate this offset from"""
        offset: IntegerSocket
        """The offset to apply to the `Index` before evaluating the input field"""

    class _Outputs(SocketAccessor):
        integer: IntegerSocket
        """The field evaluated at the offset `Index` value"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        integer: InputInteger = 0,
        index: InputInteger = 0,
        offset: InputInteger = 0,
    ):
        super().__init__(**{"Integer": integer, "Index": index, "Offset": offset})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        integer = tree.inputs.integer(
            "Integer",
            0,
            description="The field to evaluate at the given `Index` + `Offset` on the point domain",
            hide_value=True,
        )
        index = tree.inputs.integer(
            "Index",
            0,
            description="The `Index` at which to evaluate this offset from",
            min_value=0,
            default_input="INDEX",
        )
        offset = tree.inputs.integer(
            "Offset",
            0,
            description="The offset to apply to the `Index` before evaluating the input field",
        )
        integer_1 = tree.outputs.integer(
            "Integer", description="The field evaluated at the offset `Index` value"
        )

        integer.point.at(OffsetIndex(index=index, offset=offset)) >> integer_1


ASSET = OffsetInteger

ASSET_METADATA = {
    "catalog_id": "7ccb8802-a69f-483e-bf6e-4a47aaa9e940",
}
