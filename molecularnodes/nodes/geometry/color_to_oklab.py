# Node-group asset 'Color to OKLab' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    ColorSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputColor
from ._shared.oklab_matrices import OKLabMatrices


class ColorToOKLab(AssetGeometryGroup):
    """
    Color to OKLab

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
    o.oklab : VectorSocket
        OKLab
    """

    _name = "Color to OKLab"
    _asset_name = "Color to OKLab"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "COLOR"

    class _Inputs(SocketAccessor):
        color: ColorSocket
        """Color"""

    class _Outputs(SocketAccessor):
        oklab: VectorSocket
        """OKLab"""

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
        oklab = tree.outputs.vector("OKLab", subtype="XYZ")

        vector_math = g.VectorMath.power(
            g.TransformPoint(vector=color, transform=OKLabMatrices().o.m1),
            g.Value(0.3333333),
        )
        vector_math.o.vector.transform(OKLabMatrices().o.m2) >> oklab


ASSET = ColorToOKLab

ASSET_METADATA = {
    "catalog_id": "dafbb31f-8cbb-49a9-902d-470e1791a0b4",
}
