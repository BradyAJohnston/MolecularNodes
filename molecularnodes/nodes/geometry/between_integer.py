# Node-group asset "Between Integer" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputInteger


class BetweenInteger(AssetGeometryGroup):
    """
    Between Integer

    Parameters
    ----------
    value : InputInteger
        The value to test if it exists within the bounds
    lower : InputInteger
        The lower bounds for the test
    upper : InputInteger
        The upper bounds for the test

    Inputs
    ------
    i.value : IntegerSocket
        The value to test if it exists within the bounds
    i.lower : IntegerSocket
        The lower bounds for the test
    i.upper : IntegerSocket
        The upper bounds for the test

    Outputs
    -------
    o.boolean : BooleanSocket
        Whether the input `Value` is between (and including) the lower and upper bounds
    """

    _name = "Between Integer"
    _asset_name = "Between Integer"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {"node_tool_idname": "geometry.between_integer"}

    class _Inputs(SocketAccessor):
        value: IntegerSocket
        """The value to test if it exists within the bounds"""
        lower: IntegerSocket
        """The lower bounds for the test"""
        upper: IntegerSocket
        """The upper bounds for the test"""

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
        value: InputInteger = 0,
        lower: InputInteger = 0,
        upper: InputInteger = 19,
    ):
        super().__init__(**{"Value": value, "Lower": lower, "Upper": upper})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        value = tree.inputs.integer(
            "Value", 0, description="The value to test if it exists within the bounds"
        )
        lower = tree.inputs.integer(
            "Lower", 0, description="The lower bounds for the test"
        )
        upper = tree.inputs.integer(
            "Upper", 19, description="The upper bounds for the test"
        )
        boolean = tree.outputs.boolean(
            "Boolean",
            description="Whether the input `Value` is between (and including) the lower and upper bounds",
        )

        ((value >= lower) & (value <= upper)) >> boolean


ASSET = BetweenInteger

ASSET_METADATA = {
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
