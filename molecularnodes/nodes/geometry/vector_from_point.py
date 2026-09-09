# Node-group asset 'Vector from Point' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    PackageLibrary,
    RotationSocket,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputVector


class VectorFromPoint(AssetGeometryGroup):
    """
    Vector from Point

    Parameters
    ----------
    target : InputVector
        Vector that is the target
    position : InputVector
        Position of the current point

    Inputs
    ------
    i.target : VectorSocket
        Vector that is the target
    i.position : VectorSocket
        Position of the current point

    Outputs
    -------
    o.vector : VectorSocket
        Vector from the current point's position to the given vector
    o.direction : VectorSocket
        Normalized output vector
    o.length : FloatSocket
        Length of the output vector
    o.rotation : RotationSocket
        Rotation
    """

    _name = "Vector from Point"
    _asset_name = "Vector from Point"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "VECTOR"
    _tree_properties = {"node_tool_idname": "geometry.vector_from_point"}

    class _Inputs(SocketAccessor):
        target: VectorSocket
        """Vector that is the target"""
        position: VectorSocket
        """Position of the current point"""

    class _Outputs(SocketAccessor):
        vector: VectorSocket
        """Vector from the current point's position to the given vector"""
        direction: VectorSocket
        """Normalized output vector"""
        length: FloatSocket
        """Length of the output vector"""
        rotation: RotationSocket
        """Rotation"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        target: InputVector = None,
        position: InputVector = None,
    ):
        super().__init__(**{"Target": target, "Position": position})

    def _build_group(self, tree):
        target = tree.inputs.vector(
            "Target",
            (0.0, 0.0, 0.0),
            description="Vector that is the target",
            min_value=-10_000.0,
            max_value=10_000.0,
        )
        position = tree.inputs.vector(
            "Position",
            (0.0, 0.0, 0.0),
            description="Position of the current point",
            min_value=-10_000.0,
            max_value=10_000.0,
            default_input="POSITION",
        )
        vector = tree.outputs.vector(
            "Vector",
            description="Vector from the current point's position to the given vector",
        )
        direction = tree.outputs.vector(
            "Direction", description="Normalized output vector"
        )
        length = tree.outputs.float("Length", description="Length of the output vector")
        rotation = tree.outputs.rotation("Rotation")

        vector_math = target - position
        align_rotation_to_vector = g.AlignRotationToVector(vector=vector_math)
        vector_math.normalize() >> direction
        vector_math.length() >> length

        vector_math >> vector
        align_rotation_to_vector >> rotation


ASSET = VectorFromPoint

ASSET_METADATA = {
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
