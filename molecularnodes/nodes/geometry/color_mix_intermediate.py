# Node-group asset "Color Mix Intermediate" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    ColorSocket,
    FloatSocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputColor, InputFloat, InputMenu
from .color_oklab_mix import ColorOKLabMix
from .color_to_oklab import ColorToOKLab
from .oklab_to_color import OKLabToColor


class ColorMixIntermediate(AssetGeometryGroup):
    """
    Color Mix Intermediate

    Parameters
    ----------
    factor : InputFloat
        Factor
    menu : InputMenu | Literal["Linear", "OKLab"]
        Menu
    socket_2 : InputBoolean
        Intermediate
    a : InputColor
        A
    socket_4 : InputColor
        Intermediate
    b : InputColor
        B

    Inputs
    ------
    i.factor : FloatSocket
        Factor
    i.menu : MenuSocket
        Menu
    i.socket_2 : BooleanSocket
        Intermediate
    i.a : ColorSocket
        A
    i.socket_4 : ColorSocket
        Intermediate
    i.b : ColorSocket
        B

    Outputs
    -------
    o.output : ColorSocket
        Output
    """

    _name = "Color Mix Intermediate"
    _asset_name = "Color Mix Intermediate"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "COLOR"

    class _Inputs(SocketAccessor):
        factor: FloatSocket
        """Factor"""
        menu: MenuSocket
        """Menu"""
        socket_2: BooleanSocket
        """Intermediate"""
        a: ColorSocket
        """A"""
        socket_4: ColorSocket
        """Intermediate"""
        b: ColorSocket
        """B"""

    class _Outputs(SocketAccessor):
        output: ColorSocket
        """Output"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        factor: InputFloat = 0.5,
        menu: InputMenu | Literal["Linear", "OKLab"] = "Linear",
        socket_2: InputBoolean = False,
        a: InputColor = None,
        socket_4: InputColor = None,
        b: InputColor = None,
    ):
        super().__init__(
            **{"Factor": factor, "Menu": menu, "A": a, "B": b},
            _named_links=[("Intermediate", socket_2), ("Intermediate", socket_4)],
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        factor = tree.inputs.float(
            "Factor", 0.5, min_value=0.0, max_value=1.0, subtype="FACTOR"
        )
        menu = tree.inputs.menu("Menu", optional_label=True)
        intermediate = tree.inputs.boolean("Intermediate", False)
        a = tree.inputs.color("A", (0.0769495, 0.4785124, 0.5, 1.0))
        intermediate_1 = tree.inputs.color("Intermediate", (0.5, 0.5, 0.5, 1.0))
        b = tree.inputs.color("B", (0.5, 0.1594808, 0.05802507, 1.0))
        output = tree.outputs.color("Output", (0.8, 0.8, 0.8, 1.0))

        group = ColorToOKLab(color=intermediate_1)
        mix = g.Mix(
            factor_float=factor,
            a_color=a,
            b_color=intermediate_1,
            data_type="RGBA",
            clamp_factor=True,
        )
        mix_1 = g.Mix(
            factor_float=factor,
            a_color=intermediate_1,
            b_color=b,
            data_type="RGBA",
            clamp_factor=True,
        )
        mix_2 = g.Mix(
            factor_float=factor,
            a_vector=ColorToOKLab(color=a),
            b_vector=group,
            data_type="VECTOR",
            clamp_factor=True,
        )
        mix_3 = g.Mix(
            factor_float=factor,
            a_vector=group,
            b_vector=ColorToOKLab(color=b),
            data_type="VECTOR",
            clamp_factor=True,
        )
        mix_4 = g.Mix(
            factor_float=factor,
            a_vector=mix_2.o.result_vector,
            b_vector=mix_3.o.result_vector,
            data_type="VECTOR",
            clamp_factor=True,
        )
        combine_color = g.CombineColor(
            red=mix_1.o.result_color.r,
            green=mix_1.o.result_color.g,
            blue=mix_1.o.result_color.b,
            alpha=intermediate_1.a,
        )
        combine_color_1 = g.CombineColor(
            red=mix.o.result_color.r,
            green=mix.o.result_color.g,
            blue=mix.o.result_color.b,
            alpha=g.SeparateColor(color=a).o.alpha,
        )
        mix_5 = g.Mix(
            factor_float=factor,
            a_color=combine_color_1,
            b_color=combine_color,
            data_type="RGBA",
            clamp_factor=True,
        )
        combine_color_2 = g.CombineColor(
            red=mix_5.o.result_color.r,
            green=mix_5.o.result_color.g,
            blue=mix_5.o.result_color.b,
            alpha=combine_color_1.o.color.a,
        )
        result = g.Mix(
            factor_float=factor,
            a_color=a,
            b_color=b,
            data_type="RGBA",
            clamp_factor=True,
        ).o.result_color
        combine_color_3 = g.CombineColor(
            red=result.r,
            green=result.g,
            blue=result.b,
            alpha=g.SeparateColor(color=a).o.alpha,
        )
        switch = intermediate.switch.color(
            ColorOKLabMix(factor=factor, a=a, b=b),
            OKLabToColor(oklab=mix_4.o.result_vector),
        )
        (
            g.MenuSwitch.color(
                menu,
                {
                    "Linear": intermediate.switch.color(
                        combine_color_3, combine_color_2
                    ),
                    "OKLab": switch,
                },
            )
            >> output
        )

        menu.default_value = "Linear"


ASSET = ColorMixIntermediate

ASSET_METADATA = {
    "catalog_id": "d3f975df-8408-4972-a669-8187a57e01d0",
}
