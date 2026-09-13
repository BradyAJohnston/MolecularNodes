# Node-group asset "Color OKLab Mix" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    ColorSocket,
    FloatSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputColor, InputFloat
from .color_to_oklab import ColorToOKLab
from .oklab_to_color import OKLabToColor


class ColorOKLabMix(AssetGeometryGroup):
    """
    Color OKLab Mix

    Parameters
    ----------
    factor : InputFloat
        Factor
    a : InputColor
        A
    b : InputColor
        B

    Inputs
    ------
    i.factor : FloatSocket
        Factor
    i.a : ColorSocket
        A
    i.b : ColorSocket
        B

    Outputs
    -------
    o.result : ColorSocket
        Result
    """

    _name = "Color OKLab Mix"
    _asset_name = "Color OKLab Mix"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "COLOR"

    class _Inputs(SocketAccessor):
        factor: FloatSocket
        """Factor"""
        a: ColorSocket
        """A"""
        b: ColorSocket
        """B"""

    class _Outputs(SocketAccessor):
        result: ColorSocket
        """Result"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        factor: InputFloat = 0.5,
        a: InputColor = None,
        b: InputColor = None,
    ):
        super().__init__(**{"Factor": factor, "A": a, "B": b})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        factor = tree.inputs.float(
            "Factor", 0.5, min_value=0.0, max_value=1.0, subtype="FACTOR"
        )
        a = tree.inputs.color("A", (0.05120815, 0.3684095, 0.7010043, 1.0))
        b = tree.inputs.color("B", (0.6960105, 0.04253292, 0.03904183, 1.0))
        result = tree.outputs.color("Result", (0.0, 0.0, 0.0, 1.0))

        mix = g.Mix(
            factor_float=factor,
            a_vector=ColorToOKLab(color=a),
            b_vector=ColorToOKLab(color=b),
            data_type="VECTOR",
            clamp_factor=True,
        )
        OKLabToColor(oklab=mix.o.result_vector) >> result


ASSET = ColorOKLabMix

ASSET_METADATA = {
    "catalog_id": "d3f975df-8408-4972-a669-8187a57e01d0",
}
