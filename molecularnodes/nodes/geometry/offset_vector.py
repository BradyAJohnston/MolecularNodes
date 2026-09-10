# Node-group asset "Offset Vector" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
    VectorSocket,
)
from nodebpy.types import InputInteger, InputVector
from .offset_index import OffsetIndex


class OffsetVector(AssetGeometryGroup):
    """
    Offset Vector

    Parameters
    ----------
    vector : InputVector
        The field to evaluate at the given `Index` + `Offset` on the point domain
    index : InputInteger
        The `Index` at which to evaluate this offset from
    offset : InputInteger
        The offset to apply to the `Index` before evaluating the input field

    Inputs
    ------
    i.vector : VectorSocket
        The field to evaluate at the given `Index` + `Offset` on the point domain
    i.index : IntegerSocket
        The `Index` at which to evaluate this offset from
    i.offset : IntegerSocket
        The offset to apply to the `Index` before evaluating the input field

    Outputs
    -------
    o.value : VectorSocket
        The field evaluated at the offset `Index` value
    """

    _name = "Offset Vector"
    _asset_name = "Offset Vector"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {"node_tool_idname": "geometry.offset_vector"}

    class _Inputs(SocketAccessor):
        vector: VectorSocket
        """The field to evaluate at the given `Index` + `Offset` on the point domain"""
        index: IntegerSocket
        """The `Index` at which to evaluate this offset from"""
        offset: IntegerSocket
        """The offset to apply to the `Index` before evaluating the input field"""

    class _Outputs(SocketAccessor):
        value: VectorSocket
        """The field evaluated at the offset `Index` value"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        vector: InputVector = None,
        index: InputInteger = 0,
        offset: InputInteger = 0,
    ):
        super().__init__(**{"Vector": vector, "Index": index, "Offset": offset})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        vector = tree.inputs.vector(
            "Vector",
            (0.0, 0.0, 0.0),
            description="The field to evaluate at the given `Index` + `Offset` on the point domain",
            hide_value=True,
            default_input="POSITION",
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
            min_value=-2147483647,
        )
        value = tree.outputs.vector(
            "Value", description="The field evaluated at the offset `Index` value"
        )

        vector.point.at(OffsetIndex(index=index, offset=offset)) >> value


ASSET = OffsetVector

ASSET_METADATA = {
    "catalog_id": "7ccb8802-a69f-483e-bf6e-4a47aaa9e940",
}
