# Node-group asset "Color Goodsell" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    ColorSocket,
    FloatSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputColor, InputFloat
from .color_oklab_offset import ColorOKLabOffset
from .select_atomic_number import SelectAtomicNumber


class ColorGoodsell(AssetGeometryGroup):
    """
    Color Goodsell

    Parameters
    ----------
    color : InputColor
        Color to apply 'Goodsell' style colors to
    factor : InputFloat
        Amount to apply the 'Goodsell Style' coloring to
    invert : InputBoolean
        Whether to invert the darkening of the colors

    Inputs
    ------
    i.color : ColorSocket
        Color to apply 'Goodsell' style colors to
    i.factor : FloatSocket
        Amount to apply the 'Goodsell Style' coloring to
    i.invert : BooleanSocket
        Whether to invert the darkening of the colors

    Outputs
    -------
    o.color : ColorSocket
        The generated color based on the node inputs
    """

    _name = "Color Goodsell"
    _asset_name = "Color Goodsell"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "COLOR"
    _tree_properties = {"node_tool_idname": "geometry.color_goodsell"}

    class _Inputs(SocketAccessor):
        color: ColorSocket
        """Color to apply 'Goodsell' style colors to"""
        factor: FloatSocket
        """Amount to apply the 'Goodsell Style' coloring to"""
        invert: BooleanSocket
        """Whether to invert the darkening of the colors"""

    class _Outputs(SocketAccessor):
        color: ColorSocket
        """The generated color based on the node inputs"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        color: InputColor = None,
        factor: InputFloat = 0.3,
        invert: InputBoolean = False,
    ):
        super().__init__(**{"Color": color, "Factor": factor, "Invert": invert})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        color = tree.inputs.color(
            "Color",
            (0.5, 0.5, 0.5, 1.0),
            description="Color to apply 'Goodsell' style colors to",
        )
        factor = tree.inputs.float(
            "Factor",
            0.3,
            description="Amount to apply the 'Goodsell Style' coloring to",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        invert = tree.inputs.boolean(
            "Invert", False, description="Whether to invert the darkening of the colors"
        )
        color_1 = tree.outputs.color(
            "Color",
            (0.0, 0.0, 0.0, 0.0),
            description="The generated color based on the node inputs",
        )

        group = SelectAtomicNumber()
        mix = g.Mix(
            factor_float=invert.switch.boolean(
                group.o.selection, group.o.inverted
            ).switch.float(factor),
            b_float=-0.4,
            clamp_factor=True,
        )
        ColorOKLabOffset(color=color, luminance=mix.o.result_float) >> color_1


ASSET = ColorGoodsell

ASSET_METADATA = {
    "catalog_id": "d3f975df-8408-4972-a669-8187a57e01d0",
}
