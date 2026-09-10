# Node-group asset "Random Color" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
    IntegerSocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputFloat, InputInteger, InputMenu
from .lch_to_oklab import LChToOKLab
from .oklab_to_color import OKLabToColor


class RandomColor(AssetGeometryGroup):
    """
    Random Color

    Parameters
    ----------
    id : InputInteger
        ID
    color_seed : InputInteger
        Seed value for the random generation of the colors
    colorspace : InputMenu | Literal["HSL", "OKLab"]
        Colorspace
    hsl_saturation : InputFloat
        Saturlation level for the random color
    hsl_lightness : InputFloat
        Lightness value for the generated random color
    oklab_luminance : InputFloat
        OKLab Luminance
    oklab_chroma : InputFloat
        OKLab Chroma

    Inputs
    ------
    i.id : IntegerSocket
        ID
    i.color_seed : IntegerSocket
        Seed value for the random generation of the colors
    i.colorspace : MenuSocket
        Colorspace
    i.hsl_saturation : FloatSocket
        Saturlation level for the random color
    i.hsl_lightness : FloatSocket
        Lightness value for the generated random color
    i.oklab_luminance : FloatSocket
        OKLab Luminance
    i.oklab_chroma : FloatSocket
        OKLab Chroma

    Outputs
    -------
    o.color : ColorSocket
        Color
    """

    _name = "Random Color"
    _asset_name = "Random Color"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "COLOR"

    class _Inputs(SocketAccessor):
        id: IntegerSocket
        """ID"""
        color_seed: IntegerSocket
        """Seed value for the random generation of the colors"""
        colorspace: MenuSocket
        """Colorspace"""
        hsl_saturation: FloatSocket
        """Saturlation level for the random color"""
        hsl_lightness: FloatSocket
        """Lightness value for the generated random color"""
        oklab_luminance: FloatSocket
        """OKLab Luminance"""
        oklab_chroma: FloatSocket
        """OKLab Chroma"""

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
        id: InputInteger = 0,
        color_seed: InputInteger = 0,
        colorspace: InputMenu | Literal["HSL", "OKLab"] = "HSL",
        hsl_saturation: InputFloat = 0.6,
        hsl_lightness: InputFloat = 0.6,
        oklab_luminance: InputFloat = 0.9,
        oklab_chroma: InputFloat = 0.2,
    ):
        super().__init__(
            **{
                "ID": id,
                "Color Seed": color_seed,
                "Colorspace": colorspace,
                "HSL Saturation": hsl_saturation,
                "HSL Lightness": hsl_lightness,
                "OKLab Luminance": oklab_luminance,
                "OKLab Chroma": oklab_chroma,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        id = tree.inputs.integer("ID", 0, hide_value=True, default_input="ID_OR_INDEX")
        color_seed = tree.inputs.integer(
            "Color Seed",
            0,
            description="Seed value for the random generation of the colors",
            min_value=-10000,
            max_value=10000,
        )
        colorspace = tree.inputs.menu("Colorspace", expanded=True, optional_label=True)
        with tree.inputs.panel("HSL"):
            hsl_saturation = tree.inputs.float(
                "HSL Saturation",
                0.6,
                description="Saturlation level for the random color",
                min_value=0.0,
                max_value=1.0,
                subtype="FACTOR",
            )
            hsl_lightness = tree.inputs.float(
                "HSL Lightness",
                0.6,
                description="Lightness value for the generated random color",
                min_value=0.0,
                max_value=1.0,
                subtype="FACTOR",
            )
        with tree.inputs.panel("OKLab"):
            oklab_luminance = tree.inputs.float("OKLab Luminance", 0.9)
            oklab_chroma = tree.inputs.float("OKLab Chroma", 0.2)
        color = tree.outputs.color("Color", (0.8, 0.8, 0.8, 1.0))

        random_value = g.RandomValue.float(
            id=id, seed=color_seed + g.Integer(integer=0)
        )
        group = LChToOKLab(
            l=oklab_luminance,
            c=oklab_chroma,
            h=random_value.o.value.map_range(to_min=-math.pi, to_max=math.pi),
        )
        (
            g.MenuSwitch.color(
                colorspace,
                {
                    "HSL": g.CombineColor.hsl(
                        random_value, hsl_saturation, hsl_lightness
                    ),
                    "OKLab": OKLabToColor(oklab=group),
                },
            )
            >> color
        )

        colorspace.default_value = "HSL"


ASSET = RandomColor

ASSET_METADATA = {
    "catalog_id": "d3f975df-8408-4972-a669-8187a57e01d0",
}
