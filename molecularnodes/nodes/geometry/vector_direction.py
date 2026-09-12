# Node-group asset "Vector Direction" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputBoolean, InputVector


class VectorDirection(AssetGeometryGroup):
    """
    Vector Direction

    Parameters
    ----------
    normalize : InputBoolean
        Normalize
    to : InputVector
        To
    from_ : InputVector
        From

    Inputs
    ------
    i.normalize : BooleanSocket
        Normalize
    i.to : VectorSocket
        To
    i.from_ : VectorSocket
        From

    Outputs
    -------
    o.direction : VectorSocket
        Vector between the points, potentially normalized
    o.distance : FloatSocket
        Distance between the points before normalization
    """

    _name = "Vector Direction"
    _asset_name = "Vector Direction"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "VECTOR"

    class _Inputs(SocketAccessor):
        normalize: BooleanSocket
        """Normalize"""
        to: VectorSocket
        """To"""
        from_: VectorSocket
        """From"""

    class _Outputs(SocketAccessor):
        direction: VectorSocket
        """Vector between the points, potentially normalized"""
        distance: FloatSocket
        """Distance between the points before normalization"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        normalize: InputBoolean = True,
        to: InputVector = None,
        from_: InputVector = None,
    ):
        super().__init__(**{"Normalize": normalize, "To": to, "From": from_})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        normalize = tree.inputs.boolean("Normalize", True)
        to = tree.inputs.vector(
            "To", (0.0, 0.0, 0.0), min_value=-10_000.0, max_value=10_000.0
        )
        from_ = tree.inputs.vector(
            "From", (0.0, 0.0, 0.0), min_value=-10_000.0, max_value=10_000.0
        )
        direction = tree.outputs.vector(
            "Direction", description="Vector between the points, potentially normalized"
        )
        distance = tree.outputs.float(
            "Distance", description="Distance between the points before normalization"
        )

        vector_math = to - from_
        vector_math.length() >> distance
        normalize.switch.vector(vector_math, vector_math.normalize()) >> direction


ASSET = VectorDirection

ASSET_METADATA = {
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
