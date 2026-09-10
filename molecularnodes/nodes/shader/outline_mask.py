# Node-group asset 'Outline Mask' (ShaderNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy import shader as s
from nodebpy.builder import (
    AssetShaderGroup,
    FloatSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputFloat
from ._shared.edge_detection import EdgeDetection


class OutlineMask(AssetShaderGroup):
    """
    Outline Mask

    Parameters
    ----------
    threshold : InputFloat
        Threshold
    thickness : InputFloat
        Thickness

    Inputs
    ------
    i.threshold : FloatSocket
        Threshold
    i.thickness : FloatSocket
        Thickness

    Outputs
    -------
    o.outline : FloatSocket
        Outline
    """

    _name = "Outline Mask"
    _asset_name = "Outline Mask"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        threshold: FloatSocket
        """Threshold"""
        thickness: FloatSocket
        """Thickness"""

    class _Outputs(SocketAccessor):
        outline: FloatSocket
        """Outline"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        threshold: InputFloat = 0.8,
        thickness: InputFloat = 0.15,
    ):
        super().__init__(**{"Threshold": threshold, "Thickness": thickness})

    def _build_group(self, tree):
        threshold = tree.inputs.float(
            "Threshold", 0.8, min_value=0.0, max_value=10_000.0
        )
        thickness = tree.inputs.float(
            "Thickness", 0.15, min_value=0.0, max_value=10_000.0
        )
        outline = tree.outputs.float("Outline")

        camera_data = s.CameraData()
        group = EdgeDetection(offset=camera_data.o.view_distance * (thickness / 100.0))
        math_1 = g.Math.less_than(
            group.o.co_planar_delta, camera_data.o.view_distance * (threshold / 100.0)
        )
        math_2 = math_1.o.value.min(g.Math.less_than(group.o.normal_delta, 0.55)).min(
            g.Math(
                value_001=group.o.object_edge,
                value=1.0,
                operation="SUBTRACT",
                use_clamp=True,
            )
        )
        1.0 - math_2 - s.Geometry().o.backfacing >> outline


ASSET = OutlineMask

ASSET_METADATA = {
    "catalog_id": "fc8d3698-34f7-4b7e-8167-a2c0391b171b",
}
