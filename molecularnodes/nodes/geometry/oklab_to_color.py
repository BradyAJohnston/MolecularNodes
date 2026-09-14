# Node-group asset "OKLab to Color" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    ColorSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputVector
from ._shared.oklab_matrices import OKLabMatrices


class OKLabToColor(AssetGeometryGroup):
    """
    OKLab to Color

    Parameters
    ----------
    oklab : InputVector
        OKLab

    Inputs
    ------
    i.oklab : VectorSocket
        OKLab

    Outputs
    -------
    o.color : ColorSocket
        Color
    """

    _name = "OKLab to Color"
    _asset_name = "OKLab to Color"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "COLOR"

    class _Inputs(SocketAccessor):
        oklab: VectorSocket
        """OKLab"""

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
        oklab: InputVector = None,
    ):
        super().__init__(**{"OKLab": oklab})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        oklab = tree.inputs.vector("OKLab", (0.0, 0.0, 0.0))
        color = tree.outputs.color("Color", (0.0, 0.0, 0.0, 1.0))

        (
            g.VectorMath.power(
                oklab.transform(OKLabMatrices().o.m2.invert()), g.Integer(integer=3)
            ).o.vector.transform(OKLabMatrices().o.m1.invert())
            >> color
        )


ASSET = OKLabToColor

ASSET_METADATA = {
    "catalog_id": "dafbb31f-8cbb-49a9-902d-470e1791a0b4",
}
