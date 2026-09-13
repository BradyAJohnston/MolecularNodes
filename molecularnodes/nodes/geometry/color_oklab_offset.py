# Node-group asset "Color OKLab Offset" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
import math
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    ColorSocket,
    FloatSocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputColor, InputFloat, InputMenu
from .color_to_oklab import ColorToOKLab
from .oklab_offset_lch import OKLabOffsetLCh
from .oklab_to_color import OKLabToColor


class ColorOKLabOffset(AssetGeometryGroup):
    """
    Color OKLab Offset

    Parameters
    ----------
    color : InputColor
        Color
    colorspace : InputMenu | Literal["OKLab", "HSL"]
        Colorspace
    luminance : InputFloat
        Luminance
    saturation : InputFloat
        Saturation
    lightness : InputFloat
        Lightness
    hue : InputFloat
        Hue

    Inputs
    ------
    i.color : ColorSocket
        Color
    i.colorspace : MenuSocket
        Colorspace
    i.luminance : FloatSocket
        Luminance
    i.saturation : FloatSocket
        Saturation
    i.lightness : FloatSocket
        Lightness
    i.hue : FloatSocket
        Hue

    Outputs
    -------
    o.color : ColorSocket
        Color
    """

    _name = "Color OKLab Offset"
    _asset_name = "Color OKLab Offset"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "COLOR"

    class _Inputs(SocketAccessor):
        color: ColorSocket
        """Color"""
        colorspace: MenuSocket
        """Colorspace"""
        luminance: FloatSocket
        """Luminance"""
        saturation: FloatSocket
        """Saturation"""
        lightness: FloatSocket
        """Lightness"""
        hue: FloatSocket
        """Hue"""

    class _Outputs(SocketAccessor):
        color: ColorSocket
        """Color"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        color: InputColor = None,
        colorspace: InputMenu | Literal["OKLab", "HSL"] = "OKLab",
        luminance: InputFloat = 0.0,
        saturation: InputFloat = 0.0,
        lightness: InputFloat = 0.0,
        hue: InputFloat = 0.0,
    ):
        super().__init__(
            **{
                "Color": color,
                "Colorspace": colorspace,
                "Luminance": luminance,
                "Saturation": saturation,
                "Lightness": lightness,
                "Hue": hue,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        color = tree.inputs.color("Color", (0.07984344, 1.0, 0.2559341, 1.0))
        colorspace = tree.inputs.menu("Colorspace", expanded=True, optional_label=True)
        luminance = tree.inputs.float(
            "Luminance", 0.0, min_value=-1.0, max_value=1.0, subtype="FACTOR"
        )
        saturation = tree.inputs.float(
            "Saturation", 0.0, min_value=0.0, max_value=1.0, subtype="FACTOR"
        )
        lightness = tree.inputs.float(
            "Lightness", 0.0, min_value=0.0, max_value=1.0, subtype="FACTOR"
        )
        hue = tree.inputs.float(
            "Hue", 0.0, min_value=-math.pi, max_value=math.pi, subtype="FACTOR"
        )
        color_1 = tree.outputs.color("Color", (0.0, 0.0, 0.0, 1.0))

        separate_color = g.SeparateColor.hsl(color)
        combine_color = g.CombineColor.hsl(
            separate_color.o.red + hue,
            separate_color.o.green + saturation,
            separate_color.o.blue + lightness,
            separate_color.o.alpha,
        )
        (
            g.MenuSwitch.color(
                colorspace,
                {
                    "OKLab": OKLabToColor(
                        oklab=OKLabOffsetLCh(
                            oklab=ColorToOKLab(color=color), l=luminance, h=hue
                        )
                    ),
                    "HSL": combine_color,
                },
            )
            >> color_1
        )

        colorspace.default_value = "OKLab"


ASSET = ColorOKLabOffset

ASSET_METADATA = {
    "catalog_id": "d3f975df-8408-4972-a669-8187a57e01d0",
}
