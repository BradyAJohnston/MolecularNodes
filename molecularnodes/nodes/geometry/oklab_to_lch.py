# Node-group asset "OKLab to LCh" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputVector


class OKLabToLCh(AssetGeometryGroup):
    """
    OKLab to LCh

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
    o.l : FloatSocket
        L
    o.c : FloatSocket
        C
    o.h : FloatSocket
        h
    """

    _name = "OKLab to LCh"
    _asset_name = "OKLab to LCh"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"

    class _Inputs(SocketAccessor):
        oklab: VectorSocket
        """OKLab"""

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
        oklab: InputVector = None,
    ):
        super().__init__(**{"OKLab": oklab})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        oklab = tree.inputs.vector(
            "OKLab", (0.0, 0.0, 0.0), min_value=-10_000.0, max_value=10_000.0
        )
        l = tree.outputs.float("L")
        c_ = tree.outputs.float("C")
        h = tree.outputs.float("h", subtype="ANGLE")

        with g.Frame("C"):
            integer = g.Integer(integer=2)
            (oklab.y**integer + oklab.z**integer).sqrt() >> c_
        oklab.z.atan2(oklab.y) >> h

        oklab.x >> l


ASSET = OKLabToLCh

ASSET_METADATA = {
    "catalog_id": "dafbb31f-8cbb-49a9-902d-470e1791a0b4",
}
