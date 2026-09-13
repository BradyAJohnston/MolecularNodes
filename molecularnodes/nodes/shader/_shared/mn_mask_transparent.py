# Node group ".MN_mask_transparent" (ShaderNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import ShaderNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy import shader as s
from nodebpy.builder import CustomShaderGroup, FloatSocket, SocketAccessor
from nodebpy.types import InputFloat


class MN_mask_transparent(CustomShaderGroup):
    """
    .MN_mask_transparent

    Parameters
    ----------
    value : InputFloat
        Value

    Inputs
    ------
    i.value : FloatSocket
        Value

    Outputs
    -------
    o.value : FloatSocket
        Value
    """

    _name = ".MN_mask_transparent"

    class _Inputs(SocketAccessor):
        value: FloatSocket
        """Value"""

    class _Outputs(SocketAccessor):
        value: FloatSocket
        """Value"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        value: InputFloat = 0.0,
    ):
        super().__init__(**{"Value": value})

    def _build_group(self, tree: TreeBuilder[ShaderNodeTree]) -> None:
        value = tree.inputs.float("Value", 0.0, min_value=-10_000.0, max_value=10_000.0)
        value_1 = tree.outputs.float("Value")

        (
            g.Math.less_than(s.LightPath().o.transparent_depth, 1.0).o.value * value
            >> value_1
        )
