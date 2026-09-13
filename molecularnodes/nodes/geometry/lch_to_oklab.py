# Node-group asset "LCh to OKLab" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from nodebpy.types import InputFloat


class LChToOKLab(AssetGeometryGroup):
    """
    LCh to OKLab

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
    o.oklab : VectorSocket
        OKLab
    """

    _name = "LCh to OKLab"
    _asset_name = "LCh to OKLab"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"

    class _Inputs(SocketAccessor):
        l: FloatSocket
        """L"""
        c: FloatSocket
        """C"""
        h: FloatSocket

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
        l: InputFloat = 0.0,
        c: InputFloat = 0.0,
        h: InputFloat = 0.0,
    ):
        super().__init__(**{"L": l, "C": c, "h": h})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        l = tree.inputs.float("L", 0.0)
        c_ = tree.inputs.float("C", 0.0)
        h = tree.inputs.float("h", 0.0, subtype="ANGLE")
        oklab = tree.outputs.vector("OKLab")

        with g.Frame("a"):
            math_1 = c_ * h.cos()
        with g.Frame("b"):
            math_2 = c_ * h.sin()
        combine_xyz = g.CombineXYZ(x=l, y=math_1, z=math_2)

        combine_xyz >> oklab


ASSET = LChToOKLab

ASSET_METADATA = {
    "catalog_id": "dafbb31f-8cbb-49a9-902d-470e1791a0b4",
}
