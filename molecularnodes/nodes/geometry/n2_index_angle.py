# Node-group asset "2 Index Angle" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputInteger, InputVector
from .vector_angle import VectorAngle


class Group2IndexAngle(AssetGeometryGroup):
    """
    2 Index Angle

    Parameters
    ----------
    position : InputVector
        The `Position` vectors to use for the angle calculation
    index_a : InputInteger
        First end point for the angle calculation around the current point
    index_b : InputInteger
        The `Index` for the middle point in the angle calculation, defaulting to the current point
    index_c : InputInteger
        Last end point for the angle calculation around the current point

    Inputs
    ------
    i.position : VectorSocket
        The `Position` vectors to use for the angle calculation
    i.index_a : IntegerSocket
        First end point for the angle calculation around the current point
    i.index_b : IntegerSocket
        The `Index` for the middle point in the angle calculation, defaulting to the current point
    i.index_c : IntegerSocket
        Last end point for the angle calculation around the current point

    Outputs
    -------
    o.angle : FloatSocket
        Angle of the line A -> Self -> C in radians
    """

    _name = "2 Index Angle"
    _asset_name = "2 Index Angle"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {"node_tool_idname": "geometry.2_point_angle"}

    class _Inputs(SocketAccessor):
        position: VectorSocket
        """The `Position` vectors to use for the angle calculation"""
        index_a: IntegerSocket
        """First end point for the angle calculation around the current point"""
        index_b: IntegerSocket
        """The `Index` for the middle point in the angle calculation, defaulting to the current point"""
        index_c: IntegerSocket
        """Last end point for the angle calculation around the current point"""

    class _Outputs(SocketAccessor):
        angle: FloatSocket
        """Angle of the line A -> Self -> C in radians"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        position: InputVector = None,
        index_a: InputInteger = 0,
        index_b: InputInteger = 0,
        index_c: InputInteger = 2,
    ):
        super().__init__(
            **{
                "Position": position,
                "Index A": index_a,
                "Index B": index_b,
                "Index C": index_c,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        position = tree.inputs.vector(
            "Position",
            (0.0, 0.0, 0.0),
            description="The `Position` vectors to use for the angle calculation",
            default_input="POSITION",
        )
        index_a = tree.inputs.integer(
            "Index A",
            0,
            description="First end point for the angle calculation around the current point",
            min_value=0,
        )
        index_b = tree.inputs.integer(
            "Index B",
            0,
            description="The `Index` for the middle point in the angle calculation, defaulting to the current point",
            min_value=0,
            default_input="INDEX",
        )
        index_c = tree.inputs.integer(
            "Index C",
            2,
            description="Last end point for the angle calculation around the current point",
            min_value=0,
        )
        angle = tree.outputs.float(
            "Angle", description="Angle of the line A -> Self -> C in radians"
        )

        evaluate_at_index = position.point.at(index_b)
        (
            VectorAngle(
                a=position.point.at(index_a.point.at(index_b)) - evaluate_at_index,
                b=position.point.at(index_c.point.at(index_b)) - evaluate_at_index,
            )
            >> angle
        )


ASSET = Group2IndexAngle

ASSET_METADATA = {
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
