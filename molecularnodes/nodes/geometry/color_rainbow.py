# Node-group asset "Color Rainbow" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from nodebpy.types import InputFloat, InputMenu
from .chain_parameter import ChainParameter
from .lch_to_oklab import LChToOKLab
from .oklab_to_color import OKLabToColor
from .residue_parameter import ResidueParameter
from .structure_parameter import StructureParameter


class ColorRainbow(AssetGeometryGroup):
    """
    Color Rainbow

    Parameters
    ----------
    factor : InputMenu | Literal["Residue", "Chain", "Structure"]
        Factor
    color_space : InputMenu | Literal["HSV", "OKLab"]
        Color Space
    offset : InputFloat
        Offset the starting hue of the rainbow colors. HSV is 0-1 for the entire rainbow, OKLab is -Pi to Pi for the rainbow.
    hsl_saturation : InputFloat
        The `Saturation` value of the rainbow colors
    hsl_value : InputFloat
        The `Value` value of the resulting rainbow colors
    oklab_luminance : InputFloat
        OKLab Luminance
    oklab_chroma : InputFloat
        OKLab Chroma

    Inputs
    ------
    i.factor : MenuSocket
        Factor
    i.color_space : MenuSocket
        Color Space
    i.offset : FloatSocket
        Offset the starting hue of the rainbow colors. HSV is 0-1 for the entire rainbow, OKLab is -Pi to Pi for the rainbow.
    i.hsl_saturation : FloatSocket
        The `Saturation` value of the rainbow colors
    i.hsl_value : FloatSocket
        The `Value` value of the resulting rainbow colors
    i.oklab_luminance : FloatSocket
        OKLab Luminance
    i.oklab_chroma : FloatSocket
        OKLab Chroma

    Outputs
    -------
    o.color : ColorSocket
        The generated color
    """

    _name = "Color Rainbow"
    _asset_name = "Color Rainbow"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "COLOR"
    _tree_properties = {"node_tool_idname": "geometry.color_rainbow"}

    class _Inputs(SocketAccessor):
        factor: MenuSocket
        """Factor"""
        color_space: MenuSocket
        """Color Space"""
        offset: FloatSocket
        """Offset the starting hue of the rainbow colors. HSV is 0-1 for the entire rainbow, OKLab is -Pi to Pi for the rainbow."""
        hsl_saturation: FloatSocket
        """The `Saturation` value of the rainbow colors"""
        hsl_value: FloatSocket
        """The `Value` value of the resulting rainbow colors"""
        oklab_luminance: FloatSocket
        """OKLab Luminance"""
        oklab_chroma: FloatSocket
        """OKLab Chroma"""

    class _Outputs(SocketAccessor):
        color: ColorSocket
        """The generated color"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        factor: InputMenu | Literal["Residue", "Chain", "Structure"] = "Chain",
        color_space: InputMenu | Literal["HSV", "OKLab"] = "HSV",
        offset: InputFloat = 0.0,
        hsl_saturation: InputFloat = 0.8,
        hsl_value: InputFloat = 0.8,
        oklab_luminance: InputFloat = 0.94,
        oklab_chroma: InputFloat = 0.2,
    ):
        super().__init__(
            **{
                "Factor": factor,
                "Color Space": color_space,
                "Offset": offset,
                "HSL Saturation": hsl_saturation,
                "HSL Value": hsl_value,
                "OKLab Luminance": oklab_luminance,
                "OKLab Chroma": oklab_chroma,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        factor = tree.inputs.menu("Factor", optional_label=True)
        color_space = tree.inputs.menu("Color Space", optional_label=True)
        offset = tree.inputs.float(
            "Offset",
            0.0,
            description="Offset the starting hue of the rainbow colors. HSV is 0-1 for the entire rainbow, OKLab is -Pi to Pi for the rainbow.",
            min_value=-10_000.0,
            max_value=10_000.0,
        )
        with tree.inputs.panel("HSV", default_closed=True):
            hsl_saturation = tree.inputs.float(
                "HSL Saturation",
                0.8,
                description="The `Saturation` value of the rainbow colors",
                min_value=0.0,
                max_value=1.0,
                subtype="FACTOR",
            )
            hsl_value = tree.inputs.float(
                "HSL Value",
                0.8,
                description="The `Value` value of the resulting rainbow colors",
                min_value=0.0,
                max_value=1.0,
                subtype="FACTOR",
            )
        with tree.inputs.panel("OKLab", default_closed=True):
            oklab_luminance = tree.inputs.float("OKLab Luminance", 0.94)
            oklab_chroma = tree.inputs.float("OKLab Chroma", 0.2)
        color = tree.outputs.color(
            "Color", (0.8, 0.8, 0.8, 1.0), description="The generated color"
        )

        menu_switch = g.MenuSwitch.float(
            factor,
            {
                "Residue": ResidueParameter().o.factor,
                "Chain": ChainParameter().o.factor,
                "Structure": StructureParameter().o.atom_factor,
            },
        )
        combine_color = g.CombineColor.hsv(
            (offset + menu_switch.o.output).wrap(1.0, 0.0), hsl_saturation, hsl_value
        )
        group = LChToOKLab(
            l=oklab_luminance,
            c=oklab_chroma,
            h=menu_switch.o.output * math.tau + offset,
        )
        (
            g.MenuSwitch.color(
                color_space, {"HSV": combine_color, "OKLab": OKLabToColor(oklab=group)}
            )
            >> color
        )

        factor.default_value = "Chain"
        color_space.default_value = "HSV"


ASSET = ColorRainbow

ASSET_METADATA = {
    "catalog_id": "d3f975df-8408-4972-a669-8187a57e01d0",
}
