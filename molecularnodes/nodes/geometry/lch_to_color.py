# Node-group asset "LCh to Color" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    AssetGeometryGroup,
    ColorSocket,
    FloatSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputFloat
from .lch_to_oklab import LChToOKLab
from .oklab_to_color import OKLabToColor


class LChToColor(AssetGeometryGroup):
    """
    LCh to Color

    Parameters
    ----------
    l : InputFloat
        L
    c : InputFloat
        C
    h : InputFloat
        h

    Inputs
    ------
    i.l : FloatSocket
        L
    i.c : FloatSocket
        C
    i.h : FloatSocket
        h

    Outputs
    -------
    o.color : ColorSocket
        Color
    """

    _name = "LCh to Color"
    _asset_name = "LCh to Color"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "COLOR"

    class _Inputs(SocketAccessor):
        l: FloatSocket
        """L"""
        c: FloatSocket
        """C"""
        h: FloatSocket

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
        l: InputFloat = 0.0,
        c: InputFloat = 0.0,
        h: InputFloat = 0.0,
    ):
        super().__init__(**{"L": l, "C": c, "h": h})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        l = tree.inputs.float("L", 0.0)
        c_ = tree.inputs.float("C", 0.0)
        h = tree.inputs.float("h", 0.0, subtype="ANGLE")
        color = tree.outputs.color("Color", (0.0, 0.0, 0.0, 1.0))

        OKLabToColor(oklab=LChToOKLab(l=l, c=c_, h=h)) >> color


ASSET = LChToColor

ASSET_METADATA = {
    "catalog_id": "d3f975df-8408-4972-a669-8187a57e01d0",
}
