# Node-group asset 'OKLab Offset LCh' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputFloat, InputVector
from .lch_to_oklab import LChToOKLab
from .oklab_to_lch import OKLabToLCh


class OKLabOffsetLCh(AssetGeometryGroup):
    """
    OKLab Offset LCh

    Parameters
    ----------
    oklab : InputVector
        OKLab
    l : InputFloat
        L
    h : InputFloat
        h

    Inputs
    ------
    i.oklab : VectorSocket
        OKLab
    i.l : FloatSocket
        L
    i.h : FloatSocket
        h

    Outputs
    -------
    o.oklab : VectorSocket
        OKLab
    """

    _name = "OKLab Offset LCh"
    _asset_name = "OKLab Offset LCh"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"

    class _Inputs(SocketAccessor):
        oklab: VectorSocket
        """OKLab"""
        l: FloatSocket
        """L"""
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
        oklab: InputVector = None,
        l: InputFloat = 0.0,
        h: InputFloat = 0.0,
    ):
        super().__init__(**{"OKLab": oklab, "L": l, "h": h})

    def _build_group(self, tree):
        oklab = tree.inputs.vector(
            "OKLab", (0.64, 0.0, 0.0001), min_value=-10_000.0, max_value=10_000.0
        )
        l = tree.inputs.float("L", 0.0, min_value=-10_000.0, max_value=10_000.0)
        h = tree.inputs.float("h", 0.0, min_value=-10_000.0, max_value=10_000.0)
        oklab_1 = tree.outputs.vector("OKLab")

        group = OKLabToLCh(oklab=oklab)
        LChToOKLab(l=l + group.o.l, c=group.o.c, h=group.o.h + h) >> oklab_1


ASSET = OKLabOffsetLCh

ASSET_METADATA = {
    "catalog_id": "dafbb31f-8cbb-49a9-902d-470e1791a0b4",
}
