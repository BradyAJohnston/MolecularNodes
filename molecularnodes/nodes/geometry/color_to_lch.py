# Node-group asset 'Color to LCh' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy.builder import (
    AssetGeometryGroup,
    ColorSocket,
    FloatSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputColor
from .color_to_oklab import ColorToOKLab
from .oklab_to_lch import OKLabToLCh


class ColorToLCh(AssetGeometryGroup):
    """
    Color to LCh

    Parameters
    ----------
    color : InputColor
        Color

    Inputs
    ------
    i.color : ColorSocket
        Color

    Outputs
    -------
    o.l : FloatSocket
        L
    o.c : FloatSocket
        C
    o.h : FloatSocket
        h
    """

    _name = "Color to LCh"
    _asset_name = "Color to LCh"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "COLOR"

    class _Inputs(SocketAccessor):
        color: ColorSocket
        """Color"""

    class _Outputs(SocketAccessor):
        l: FloatSocket
        """L"""
        c: FloatSocket
        """C"""
        h: FloatSocket

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        color: InputColor = None,
    ):
        super().__init__(**{"Color": color})

    def _build_group(self, tree):
        color = tree.inputs.color("Color", (0.0, 0.0, 0.0, 1.0))
        l = tree.outputs.float("L")
        c_ = tree.outputs.float("C")
        h = tree.outputs.float("h", subtype="ANGLE")

        group = OKLabToLCh(oklab=ColorToOKLab(color=color))

        group >> l
        group.o.c >> c_
        group.o.h >> h


ASSET = ColorToLCh

ASSET_METADATA = {
    "catalog_id": "d3f975df-8408-4972-a669-8187a57e01d0",
}
