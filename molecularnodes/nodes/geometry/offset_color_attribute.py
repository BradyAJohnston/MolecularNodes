# Node-group asset "Offset Color Attribute" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    AssetGeometryGroup,
    ColorSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputInteger
from .color import Color
from .offset_index import OffsetIndex


class OffsetColorAttribute(AssetGeometryGroup):
    """
    Offset Color Attribute

    Parameters
    ----------
    index : InputInteger
        The `Index` at which to evaluate this offset from
    offset : InputInteger
        The offset to apply to the `Index` before evaluating the input field

    Inputs
    ------
    i.index : IntegerSocket
        The `Index` at which to evaluate this offset from
    i.offset : IntegerSocket
        The offset to apply to the `Index` before evaluating the input field

    Outputs
    -------
    o.color : ColorSocket
        The field evaluated at the offset `Index` value
    """

    _name = "Offset Color Attribute"
    _asset_name = "Offset Color Attribute"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {"node_tool_idname": "geometry.offset_color_attribute"}

    class _Inputs(SocketAccessor):
        index: IntegerSocket
        """The `Index` at which to evaluate this offset from"""
        offset: IntegerSocket
        """The offset to apply to the `Index` before evaluating the input field"""

    class _Outputs(SocketAccessor):
        color: ColorSocket
        """The field evaluated at the offset `Index` value"""

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

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        index = tree.inputs.integer(
            "Index",
            0,
            description="The `Index` at which to evaluate this offset from",
            default_input="INDEX",
        )
        offset = tree.inputs.integer(
            "Offset",
            0,
            description="The offset to apply to the `Index` before evaluating the input field",
        )
        color = tree.outputs.color(
            "Color",
            (0.0, 0.0, 0.0, 0.0),
            description="The field evaluated at the offset `Index` value",
        )

        Color(index=OffsetIndex(index=index, offset=offset)) >> color


ASSET = OffsetColorAttribute

ASSET_METADATA = {
    "catalog_id": "7ccb8802-a69f-483e-bf6e-4a47aaa9e940",
}
