# Node-group asset "Force Brownian" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from nodebpy.types import InputFloat, InputVector


class ForceBrownian(AssetGeometryGroup):
    """
    Force Brownian

    Parameters
    ----------
    add : InputVector
        Add
    small : InputFloat
        Small
    large : InputFloat
        Large

    Inputs
    ------
    i.add : VectorSocket
        Add
    i.small : FloatSocket
        Small
    i.large : FloatSocket
        Large

    Outputs
    -------
    o.force : VectorSocket
        Force
    """

    _name = "Force Brownian"
    _asset_name = "Force Brownian"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        add: VectorSocket
        """Add"""
        small: FloatSocket
        """Small"""
        large: FloatSocket
        """Large"""

    class _Outputs(SocketAccessor):
        force: VectorSocket
        """Force"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        add: InputVector = None,
        small: InputFloat = 0.1,
        large: InputFloat = 0.0,
    ):
        super().__init__(**{"Add": add, "Small": small, "Large": large})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        add = tree.inputs.vector(
            "Add",
            (0.0, 0.0, 0.0),
            min_value=-10_000.0,
            max_value=10_000.0,
            hide_value=True,
        )
        small = tree.inputs.float("Small", 0.1, min_value=0.0, max_value=10_000.0)
        large = tree.inputs.float("Large", 0.0, min_value=0.0, max_value=10_000.0)
        force = tree.outputs.vector("Force")

        scene_time = g.SceneTime()
        noise_texture = g.NoiseTexture(
            w=scene_time.o.seconds,
            scale=20.0,
            detail=15.0,
            roughness=1.0,
            noise_dimensions="4D",
        )
        noise_texture_1 = g.NoiseTexture(
            w=scene_time.o.seconds,
            scale=0.1,
            detail=15.0,
            roughness=1.0,
            distortion=10.0,
            noise_dimensions="4D",
        )
        (
            g.VectorMath.scale(noise_texture.o.color, small / 10.0).o.vector
            + (add + g.VectorMath.scale(noise_texture_1.o.color, large / 50.0).o.vector)
            >> force
        )


ASSET = ForceBrownian

ASSET_METADATA = {
    "catalog_id": "c2c958af-5095-4fc2-884d-709bba965fc4",
}
