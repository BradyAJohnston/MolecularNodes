# Node-group asset 'Index Distance' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    IntegerSocket,
    PackageLibrary,
    RotationSocket,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputInteger, InputVector
from .vector_from_point import VectorFromPoint


class IndexDistance(AssetGeometryGroup):
    """
    Index Distance

    Parameters
    ----------
    index : InputInteger
        Index
    target_index : InputInteger
        Index for the selected point to measure to
    position : InputVector
        Position

    Inputs
    ------
    i.index : IntegerSocket
        Index
    i.target_index : IntegerSocket
        Index for the selected point to measure to
    i.position : VectorSocket
        Position

    Outputs
    -------
    o.vector : VectorSocket
        Vector from the current point to the indexed point
    o.direction : VectorSocket
        Normalized vector from the current point to the indexed point
    o.distance : FloatSocket
        Distance from the current point to the indexed point
    o.rotation : RotationSocket
        Rotation
    """

    _name = "Index Distance"
    _asset_name = "Index Distance"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {"node_tool_idname": "geometry.point_distance"}

    class _Inputs(SocketAccessor):
        index: IntegerSocket
        """Index"""
        target_index: IntegerSocket
        """Index for the selected point to measure to"""
        position: VectorSocket
        """Position"""

    class _Outputs(SocketAccessor):
        vector: VectorSocket
        """Vector from the current point to the indexed point"""
        direction: VectorSocket
        """Normalized vector from the current point to the indexed point"""
        distance: FloatSocket
        """Distance from the current point to the indexed point"""
        rotation: RotationSocket
        """Rotation"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        index: InputInteger = 0,
        target_index: InputInteger = 100,
        position: InputVector = None,
    ):
        super().__init__(
            **{"Index": index, "Target Index": target_index, "Position": position}
        )

    def _build_group(self, tree):
        index = tree.inputs.integer("Index", 0, min_value=0, default_input="INDEX")
        target_index = tree.inputs.integer(
            "Target Index",
            100,
            description="Index for the selected point to measure to",
            min_value=0,
        )
        position = tree.inputs.vector(
            "Position", (0.0, 0.0, 0.0), default_input="POSITION"
        )
        vector = tree.outputs.vector(
            "Vector", description="Vector from the current point to the indexed point"
        )
        direction = tree.outputs.vector(
            "Direction",
            description="Normalized vector from the current point to the indexed point",
        )
        distance = tree.outputs.float(
            "Distance",
            description="Distance from the current point to the indexed point",
        )
        rotation = tree.outputs.rotation("Rotation")

        group = VectorFromPoint(
            target=position.point.at(target_index.point.at(index)),
            position=g.Position().o.position.point.at(index),
        )

        group >> vector
        group.o.direction >> direction
        group.o.length >> distance
        group.o.rotation >> rotation


ASSET = IndexDistance

ASSET_METADATA = {
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
