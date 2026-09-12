# Node-group asset "Between Float" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputFloat


class BetweenFloat(AssetGeometryGroup):
    """
    Between Float

    Parameters
    ----------
    value : InputFloat
        The value to test if it is between the lower and upper bounds
    lower : InputFloat
        Test if the `Value` is greater than or equal to this
    upper : InputFloat
        Test if the `Value` is less than or equal to this

    Inputs
    ------
    i.value : FloatSocket
        The value to test if it is between the lower and upper bounds
    i.lower : FloatSocket
        Test if the `Value` is greater than or equal to this
    i.upper : FloatSocket
        Test if the `Value` is less than or equal to this

    Outputs
    -------
    o.boolean : BooleanSocket
        Whether the input `Value` is between (and including) the lower and upper bounds
    """

    _name = "Between Float"
    _asset_name = "Between Float"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {"node_tool_idname": "geometry.between_float"}

    class _Inputs(SocketAccessor):
        value: FloatSocket
        """The value to test if it is between the lower and upper bounds"""
        lower: FloatSocket
        """Test if the `Value` is greater than or equal to this"""
        upper: FloatSocket
        """Test if the `Value` is less than or equal to this"""

    class _Outputs(SocketAccessor):
        boolean: BooleanSocket
        """Whether the input `Value` is between (and including) the lower and upper bounds"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        value: InputFloat = 0.0,
        lower: InputFloat = 0.0,
        upper: InputFloat = 0.0,
    ):
        super().__init__(**{"Value": value, "Lower": lower, "Upper": upper})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        value = tree.inputs.float(
            "Value",
            0.0,
            description="The value to test if it is between the lower and upper bounds",
        )
        lower = tree.inputs.float(
            "Lower",
            0.0,
            description="Test if the `Value` is greater than or equal to this",
        )
        upper = tree.inputs.float(
            "Upper",
            0.0,
            description="Test if the `Value` is less than or equal to this",
        )
        boolean = tree.outputs.boolean(
            "Boolean",
            description="Whether the input `Value` is between (and including) the lower and upper bounds",
        )

        ((value >= lower) & (value <= upper)) >> boolean


ASSET = BetweenFloat

ASSET_METADATA = {
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
