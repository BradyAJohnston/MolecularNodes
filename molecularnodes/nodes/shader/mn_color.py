# Node-group asset "MN Color" (ShaderNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import ShaderNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy import shader as s
from nodebpy.builder import (
    AssetShaderGroup,
    ColorSocket,
    FloatSocket,
    PackageLibrary,
    SocketAccessor,
)


class MNColor(AssetShaderGroup):
    """
    MN Color

    Outputs
    -------
    o.color : ColorSocket
        Color
    o.alpha : FloatSocket
        Alpha
    """

    _name = "MN Color"
    _asset_name = "MN Color"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "COLOR"

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        color: ColorSocket
        """Color"""
        alpha: FloatSocket
        """Alpha"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self):
        super().__init__()

    def _build_group(self, tree: TreeBuilder[ShaderNodeTree]) -> None:
        color = tree.outputs.color("Color", (0.0, 0.0, 0.0, 0.0))
        alpha = tree.outputs.float("Alpha")

        attribute = s.Attribute(attribute_name="Color")
        attribute_1 = s.Attribute(attribute_type="INSTANCER", attribute_name="Color")
        attribute_2 = s.Attribute(
            attribute_type="INSTANCER", attribute_name="is_instanced"
        )
        mix = g.Mix(
            factor_float=attribute_2.o.factor,
            a_color=attribute.o.color,
            b_color=attribute_1.o.color,
            data_type="RGBA",
            clamp_factor=True,
        )
        mix_1 = g.Mix(
            factor_float=attribute_2.o.factor,
            a_float=attribute.o.alpha,
            b_float=attribute_1.o.alpha,
            clamp_factor=True,
        )

        mix.o.result_color >> color
        mix_1 >> alpha


ASSET = MNColor

ASSET_METADATA = {
    "catalog_id": "d3f975df-8408-4972-a669-8187a57e01d0",
}
