# Node-group asset "Color Attribute Random" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
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
    StringSocket,
)
from nodebpy.types import InputFloat, InputInteger, InputMenu, InputString
from .random_color import RandomColor


class ColorAttributeRandom(AssetGeometryGroup):
    """
    Color Attribute Random

    Parameters
    ----------
    name : InputString
        Attribute to base the random color generation on
    colorspace : InputMenu | Literal["HSL", "OKLab"]
        Colorspace
    color_seed : InputInteger
        Seed value for the random generation of the colors
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
    i.name : StringSocket
        Attribute to base the random color generation on
    i.colorspace : MenuSocket
        Colorspace
    i.color_seed : IntegerSocket
        Seed value for the random generation of the colors
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
        The randomly generated color based on the input attribute
    """

    _name = "Color Attribute Random"
    _asset_name = "Color Attribute Random"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "COLOR"
    _tree_properties = {"node_tool_idname": "geometry.color_attribute_random"}

    class _Inputs(SocketAccessor):
        name: StringSocket
        """Attribute to base the random color generation on"""
        colorspace: MenuSocket
        """Colorspace"""
        color_seed: IntegerSocket
        """Seed value for the random generation of the colors"""
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
        """The randomly generated color based on the input attribute"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        name: InputString = "chain_id",
        colorspace: InputMenu | Literal["HSL", "OKLab"] = "HSL",
        color_seed: InputInteger = 0,
        hsl_saturation: InputFloat = 0.6,
        hsl_lightness: InputFloat = 0.6,
        oklab_luminance: InputFloat = 0.9,
        oklab_chroma: InputFloat = 0.2,
    ):
        super().__init__(
            **{
                "Name": name,
                "Colorspace": colorspace,
                "Color Seed": color_seed,
                "HSL Saturation": hsl_saturation,
                "HSL Lightness": hsl_lightness,
                "OKLab Luminance": oklab_luminance,
                "OKLab Chroma": oklab_chroma,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        name = tree.inputs.string(
            "Name",
            "chain_id",
            description="Attribute to base the random color generation on ",
            optional_label=True,
        )
        colorspace = tree.inputs.menu("Colorspace", optional_label=True)
        color_seed = tree.inputs.integer(
            "Color Seed",
            0,
            description="Seed value for the random generation of the colors",
            min_value=-10000,
            max_value=10000,
        )
        with tree.inputs.panel("HSL", default_closed=True):
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
        color = tree.outputs.color(
            "Color",
            (0.0, 0.0, 0.0, 0.0),
            description="The randomly generated color based on the input attribute",
        )

        _group = RandomColor(colorspace="OKLab")
        (
            RandomColor(
                id=g.NamedAttribute.integer(name).o.attribute,
                color_seed=color_seed,
                colorspace=colorspace,
                hsl_saturation=hsl_saturation,
                hsl_lightness=hsl_lightness,
                oklab_luminance=oklab_luminance,
                oklab_chroma=oklab_chroma,
            )
            >> color
        )

        colorspace.default_value = "HSL"


ASSET = ColorAttributeRandom

ASSET_METADATA = {
    "catalog_id": "d3f975df-8408-4972-a669-8187a57e01d0",
}
