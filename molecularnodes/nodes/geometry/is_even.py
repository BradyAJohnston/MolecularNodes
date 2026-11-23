# Node-group asset "Is Even" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import Default, InputInteger


class IsEven(AssetGeometryGroup):
    """
    Is Even

    Parameters
    ----------
    value : InputInteger
        Value. When unconnected: The index from the context.

    Inputs
    ------
    i.value : IntegerSocket
        Value

    Outputs
    -------
    o.even : BooleanSocket
        Even
    o.odd : BooleanSocket
        Odd
    """

    _name = "Is Even"
    _asset_name = "Is Even"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"

    class _Inputs(SocketAccessor):
        value: IntegerSocket
        """Value"""

    class _Outputs(SocketAccessor):
        even: BooleanSocket
        """Even"""
        odd: BooleanSocket
        """Odd"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        value: InputInteger = Default.INDEX,
    ):
        super().__init__(**{"Value": value})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        value = tree.inputs.integer("Value", 0, default_input="INDEX")
        even = tree.outputs.boolean("Even")
        odd = tree.outputs.boolean("Odd")

        integer_math = g.IntegerMath.modulo(value, 2)
        boolean_math = g.BooleanMath.l_not(integer_math)

        boolean_math >> even
        integer_math >> odd


ASSET = IsEven

ASSET_METADATA = {
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
