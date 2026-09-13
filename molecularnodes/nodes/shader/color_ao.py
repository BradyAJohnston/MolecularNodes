# Node-group asset "Color AO" (ShaderNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import ShaderNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy import shader as s
from nodebpy.builder import (
    AssetShaderGroup,
    ColorSocket,
    FloatSocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputColor, InputFloat, InputMenu


class ColorAO(AssetShaderGroup):
    """
    Color AO

    Parameters
    ----------
    color : InputColor
        Color
    menu : InputMenu | Literal["AO", "None"]
        Menu
    distance : InputFloat
        Distance
    exponent : InputFloat
        Exponent

    Inputs
    ------
    i.color : ColorSocket
        Color
    i.menu : MenuSocket
        Menu
    i.distance : FloatSocket
        Distance
    i.exponent : FloatSocket
        Exponent

    Outputs
    -------
    o.result : ColorSocket
        Result
    """

    _name = "Color AO"
    _asset_name = "Color AO"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "COLOR"

    class _Inputs(SocketAccessor):
        color: ColorSocket
        """Color"""
        menu: MenuSocket
        """Menu"""
        distance: FloatSocket
        """Distance"""
        exponent: FloatSocket
        """Exponent"""

    class _Outputs(SocketAccessor):
        result: ColorSocket
        """Result"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        color: InputColor = None,
        menu: InputMenu | Literal["AO", "None"] = "AO",
        distance: InputFloat = 1.0,
        exponent: InputFloat = 2.0,
    ):
        super().__init__(
            **{"Color": color, "Menu": menu, "Distance": distance, "Exponent": exponent}
        )

    def _build_group(self, tree: TreeBuilder[ShaderNodeTree]) -> None:
        color = tree.inputs.color("Color", (0.0, 0.0, 0.0, 0.0))
        menu = tree.inputs.menu("Menu", expanded=True, optional_label=True)
        distance = tree.inputs.float("Distance", 1.0, min_value=0.0, max_value=1000.0)
        exponent = tree.inputs.float("Exponent", 2.0, min_value=0.0, max_value=10_000.0)
        result = tree.outputs.color("Result", (0.8, 0.8, 0.8, 1.0))

        math_1 = g.Math(
            value_001=s.AmbientOcclusion(
                color=color, distance=distance, samples=16
            ).o.ao
            ** exponent,
            value=1.0,
            operation="SUBTRACT",
            use_clamp=True,
        )
        mix = g.Mix(
            factor_float=s.MenuSwitch.float(menu, {"AO": math_1, "None": 0.0}).o.output,
            a_color=color,
            b_color=(0.0, 0.0, 0.0, 1.0),
            data_type="RGBA",
            blend_type="MULTIPLY",
            clamp_factor=True,
        )

        mix.o.result_color >> result

        menu.default_value = "AO"


ASSET = ColorAO
