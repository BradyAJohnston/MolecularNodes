# Node-group asset "Color Backbone" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    AssetGeometryGroup,
    ColorSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputColor
from .color import Color
from .is_alpha_carbon import IsAlphaCarbon
from .is_backbone import IsBackbone
from .is_side_chain import IsSideChain


class ColorBackbone(AssetGeometryGroup):
    """
    Color Backbone

    Parameters
    ----------
    backbone : InputColor
        Backbone
    side_chain : InputColor
        Side Chain

    Inputs
    ------
    i.backbone : ColorSocket
        Backbone
    i.side_chain : ColorSocket
        Side Chain

    Outputs
    -------
    o.color : ColorSocket
        Color
    """

    _name = "Color Backbone"
    _asset_name = "Color Backbone"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "COLOR"
    _tree_properties = {"node_tool_idname": "geometry.color_backbone"}

    class _Inputs(SocketAccessor):
        backbone: ColorSocket
        """Backbone"""
        side_chain: ColorSocket
        """Side Chain"""

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
        backbone: InputColor = None,
        side_chain: InputColor = None,
    ):
        super().__init__(**{"Backbone": backbone, "Side Chain": side_chain})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        backbone = tree.inputs.color("Backbone", (0.469481, 0.24, 0.6, 1.0))
        side_chain = tree.inputs.color("Side Chain", (0.525519, 0.6, 0.24, 1.0))
        color = tree.outputs.color("Color", (0.0, 0.0, 0.0, 0.0))

        (
            (~IsAlphaCarbon().o.selection & IsSideChain().o.selection).switch.color(
                IsBackbone().o.selection.switch.color(Color(), backbone), side_chain
            )
            >> color
        )


ASSET = ColorBackbone

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
