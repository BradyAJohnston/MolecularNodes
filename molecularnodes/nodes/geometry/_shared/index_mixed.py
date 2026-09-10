# Node group 'Index Mixed' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    CustomGeometryGroup,
    FloatSocket,
    IntegerSocket,
    SocketAccessor,
)
from nodebpy.types import InputFloat, InputInteger


class IndexMixed(CustomGeometryGroup):
    """
    Index Mixed

    Parameters
    ----------
    index : InputInteger
        The `Index` at which to add the `Offset` value to
    offset : InputFloat
        The offset value to add to the to the `Index`

    Inputs
    ------
    i.index : IntegerSocket
        The `Index` at which to add the `Offset` value to
    i.offset : FloatSocket
        The offset value to add to the to the `Index`

    Outputs
    -------
    o.mixed : FloatSocket
        The sum of the `Offset` and the `Index`
    o.floor : IntegerSocket
        The floor of the `Mixed` output
    o.ceiling : IntegerSocket
        The `Ceiling` of the mixed output
    """

    _name = "Index Mixed"
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.index_mixed"}

    class _Inputs(SocketAccessor):
        index: IntegerSocket
        """The `Index` at which to add the `Offset` value to"""
        offset: FloatSocket
        """The offset value to add to the to the `Index`"""

    class _Outputs(SocketAccessor):
        mixed: FloatSocket
        """The sum of the `Offset` and the `Index`"""
        floor: IntegerSocket
        """The floor of the `Mixed` output"""
        ceiling: IntegerSocket
        """The `Ceiling` of the mixed output"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        index: InputInteger = 0,
        offset: InputFloat = 0.5,
    ):
        super().__init__(**{"Index": index, "Offset": offset})

    def _build_group(self, tree):
        index = tree.inputs.integer(
            "Index",
            0,
            description="The `Index` at which to add the `Offset` value to",
            default_input="INDEX",
        )
        offset = tree.inputs.float(
            "Offset",
            0.5,
            description="The offset value to add to the to the `Index`",
            min_value=-10_000.0,
            max_value=10_000.0,
        )
        mixed = tree.outputs.float(
            "Mixed", description="The sum of the `Offset` and the `Index`"
        )
        floor = tree.outputs.integer(
            "Floor", description="The floor of the `Mixed` output"
        )
        ceiling = tree.outputs.integer(
            "Ceiling", description="The `Ceiling` of the mixed output"
        )

        math_1 = index + offset
        float_to_integer = g.FloatToInteger(float=math_1, rounding_mode="FLOOR")
        float_to_integer.o.integer + 1.0 >> ceiling

        math_1 >> mixed
        float_to_integer >> floor
