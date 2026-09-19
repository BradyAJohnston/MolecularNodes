# Node group "Ease Switch" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import CustomGeometryGroup, FloatSocket, MenuSocket, SocketAccessor
from nodebpy.types import InputFloat, InputMenu


class EaseSwitch(CustomGeometryGroup):
    """
    Ease Switch

    Parameters
    ----------
    ease : InputMenu | Literal["In", "Out", "In Out"]
        Apply the curve at the start (In), the end (Out) or both ends (In Out) of the transition
    in_ : InputFloat
        Value used when easing In
    out : InputFloat
        Value used when easing Out
    in_out : InputFloat
        Value used when easing In Out

    Inputs
    ------
    i.ease : MenuSocket
        Apply the curve at the start (In), the end (Out) or both ends (In Out) of the transition
    i.in_ : FloatSocket
        Value used when easing In
    i.out : FloatSocket
        Value used when easing Out
    i.in_out : FloatSocket
        Value used when easing In Out

    Outputs
    -------
    o.value : FloatSocket
        Value
    """

    _name = "Ease Switch"
    _color_tag = "CONVERTER"

    class _Inputs(SocketAccessor):
        ease: MenuSocket
        """Apply the curve at the start (In), the end (Out) or both ends (In Out) of the transition"""
        in_: FloatSocket
        """Value used when easing In"""
        out: FloatSocket
        """Value used when easing Out"""
        in_out: FloatSocket
        """Value used when easing In Out"""

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
        ease: InputMenu | Literal["In", "Out", "In Out"] = "In Out",
        in_: InputFloat = 0.0,
        out: InputFloat = 0.0,
        in_out: InputFloat = 0.0,
    ):
        super().__init__(**{"Ease": ease, "In": in_, "Out": out, "In Out": in_out})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        ease = tree.inputs.menu(
            "Ease",
            description="Apply the curve at the start (In), the end (Out) or both ends (In Out) of the transition",
        )
        in_ = tree.inputs.float("In", 0.0, description="Value used when easing In")
        out = tree.inputs.float("Out", 0.0, description="Value used when easing Out")
        in_out = tree.inputs.float(
            "In Out", 0.0, description="Value used when easing In Out"
        )
        value = tree.outputs.float("Value")

        g.MenuSwitch.float(ease, {"In": in_, "Out": out, "In Out": in_out}) >> value

        ease.default_value = "In Out"
