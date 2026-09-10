# Node-group asset "Between Vector" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputVector


class BetweenVector(AssetGeometryGroup):
    """
    Between Vector

    Parameters
    ----------
    value : InputVector
        The value to test element-wise
    lower : InputVector
        The lower bounds (including) for the comparison
    upper : InputVector
        The upper bounds (including) for the comparison

    Inputs
    ------
    i.value : VectorSocket
        The value to test element-wise
    i.lower : VectorSocket
        The lower bounds (including) for the comparison
    i.upper : VectorSocket
        The upper bounds (including) for the comparison

    Outputs
    -------
    o.boolean : BooleanSocket
        If the value is between (and including) the lower the upper bounds
    """

    _name = "Between Vector"
    _asset_name = "Between Vector"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {"node_tool_idname": "geometry.between_vector"}

    class _Inputs(SocketAccessor):
        value: VectorSocket
        """The value to test element-wise"""
        lower: VectorSocket
        """The lower bounds (including) for the comparison"""
        upper: VectorSocket
        """The upper bounds (including) for the comparison"""

    class _Outputs(SocketAccessor):
        boolean: BooleanSocket
        """If the value is between (and including) the lower the upper bounds"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        value: InputVector = None,
        lower: InputVector = None,
        upper: InputVector = None,
    ):
        super().__init__(**{"Value": value, "Lower": lower, "Upper": upper})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        value = tree.inputs.vector(
            "Value", (0.0, 0.0, 0.0), description="The value to test element-wise"
        )
        lower = tree.inputs.vector(
            "Lower",
            (0.0, 0.0, 0.0),
            description="The lower bounds (including) for the comparison",
        )
        upper = tree.inputs.vector(
            "Upper",
            (0.0, 0.0, 0.0),
            description="The upper bounds (including) for the comparison",
        )
        boolean = tree.outputs.boolean(
            "Boolean",
            description="If the value is between (and including) the lower the upper bounds",
        )

        ((value >= lower) & (value <= upper)) >> boolean


ASSET = BetweenVector

ASSET_METADATA = {
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
