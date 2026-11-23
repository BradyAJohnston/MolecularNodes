# Node-group asset "Decay" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputFloat


class Decay(AssetGeometryGroup):
    """
    Decay

    Parameters
    ----------
    delta_time : InputFloat
        Delta Time
    half_life : InputFloat
        Seconds for the remaining gap to the best view to halve

    Inputs
    ------
    i.delta_time : FloatSocket
        Delta Time
    i.half_life : FloatSocket
        Seconds for the remaining gap to the best view to halve

    Outputs
    -------
    o.decay : FloatSocket
        Decay
    """

    _name = "Decay"
    _asset_name = "Decay"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"

    class _Inputs(SocketAccessor):
        delta_time: FloatSocket
        """Delta Time"""
        half_life: FloatSocket
        """Seconds for the remaining gap to the best view to halve"""

    class _Outputs(SocketAccessor):
        decay: FloatSocket
        """Decay"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        delta_time: InputFloat = 0.0,
        half_life: InputFloat = 0.25,
    ):
        super().__init__(**{"Delta Time": delta_time, "Half Life": half_life})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        delta_time = tree.inputs.float("Delta Time", 0.0)
        half_life = tree.inputs.float(
            "Half Life",
            0.25,
            description="Seconds for the remaining gap to the best view to halve",
            min_value=0.0,
            subtype="TIME_ABSOLUTE",
        )
        decay = tree.outputs.float("Decay")

        1.0 - 2.0 ** (delta_time * -1.0 / half_life.max(0.000001)) >> decay


ASSET = Decay

ASSET_METADATA = {
    "catalog_id": "85730213-4c2e-469f-b333-52ac53adf274",
}
