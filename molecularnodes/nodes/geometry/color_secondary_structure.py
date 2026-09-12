# Node-group asset "Color Secondary Structure" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from .is_helix import IsHelix
from .is_loop import IsLoop
from .is_sheet import IsSheet


class ColorSecondaryStructure(AssetGeometryGroup):
    """
    Color Secondary Structure

    Parameters
    ----------
    helix : InputColor
        Color to set for alpha helices
    sheet : InputColor
        Color to set for beta-sheets
    loop : InputColor
        Color to set for loops in the structure
    other : InputColor
        Other

    Inputs
    ------
    i.helix : ColorSocket
        Color to set for alpha helices
    i.sheet : ColorSocket
        Color to set for beta-sheets
    i.loop : ColorSocket
        Color to set for loops in the structure
    i.other : ColorSocket
        Other

    Outputs
    -------
    o.color : ColorSocket
        The colors based on secondary structure
    """

    _name = "Color Secondary Structure"
    _asset_name = "Color Secondary Structure"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "COLOR"
    _tree_properties = {"node_tool_idname": "geometry.color_sec_struct"}

    class _Inputs(SocketAccessor):
        helix: ColorSocket
        """Color to set for alpha helices"""
        sheet: ColorSocket
        """Color to set for beta-sheets"""
        loop: ColorSocket
        """Color to set for loops in the structure"""
        other: ColorSocket
        """Other"""

    class _Outputs(SocketAccessor):
        color: ColorSocket
        """The colors based on secondary structure"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        helix: InputColor = None,
        sheet: InputColor = None,
        loop: InputColor = None,
        other: InputColor = None,
    ):
        super().__init__(
            **{"Helix": helix, "Sheet": sheet, "Loop": loop, "Other": other}
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        helix = tree.inputs.color(
            "Helix",
            (0.16202937, 0.6239606, 0.19461784, 1.0),
            description="Color to set for alpha helices",
        )
        sheet = tree.inputs.color(
            "Sheet",
            (1.0, 0.1499598, 0.1499598, 1.0),
            description="Color to set for beta-sheets",
        )
        loop = tree.inputs.color(
            "Loop",
            (0.17144114, 0.3662526, 0.799103, 1.0),
            description="Color to set for loops in the structure",
        )
        other = tree.inputs.color("Other", (0.8, 0.8, 0.8, 1.0))
        color = tree.outputs.color(
            "Color",
            (0.0, 0.0, 0.0, 0.0),
            description="The colors based on secondary structure",
        )

        switch = IsSheet().o.selection.switch.color(
            IsHelix().o.selection.switch.color(other, helix), sheet
        )
        IsLoop().o.selection.switch.color(switch, loop) >> color


ASSET = ColorSecondaryStructure

ASSET_METADATA = {
    "catalog_id": "d3f975df-8408-4972-a669-8187a57e01d0",
}
