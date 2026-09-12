# Node-group asset "Vector Angle" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputVector


class VectorAngle(AssetGeometryGroup):
    """
    The angle between two vectors, in radians

    Parameters
    ----------
    a : InputVector
        The first vector for angle calculation
    b : InputVector
        The second vector for the angle calculation

    Inputs
    ------
    i.a : VectorSocket
        The first vector for angle calculation
    i.b : VectorSocket
        The second vector for the angle calculation

    Outputs
    -------
    o.angle : FloatSocket
        The angle between the two vectors in radians
    o.a_b : VectorSocket
        Axis around which the angle rotates (cross product of A and B)
    """

    _name = "Vector Angle"
    _asset_name = "Vector Angle"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "VECTOR"
    _tree_properties = {
        "description": "The angle between two vectors, in radians",
        "node_tool_idname": "geometry.vector_angle",
    }

    class _Inputs(SocketAccessor):
        a: VectorSocket
        """The first vector for angle calculation"""
        b: VectorSocket
        """The second vector for the angle calculation"""

    class _Outputs(SocketAccessor):
        angle: FloatSocket
        """The angle between the two vectors in radians"""
        a_b: VectorSocket
        """Axis around which the angle rotates (cross product of A and B)"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        a: InputVector = None,
        b: InputVector = None,
    ):
        super().__init__(**{"A": a, "B": b})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        a = tree.inputs.vector(
            "A",
            (0.0, 0.0, 0.0),
            description="The first vector for angle calculation",
            min_value=-10_000.0,
            max_value=10_000.0,
        )
        b = tree.inputs.vector(
            "B",
            (0.0, 0.0, 0.0),
            description="The second vector for the angle calculation",
            min_value=-10_000.0,
            max_value=10_000.0,
        )
        angle = tree.outputs.float(
            "Angle",
            description="The angle between the two vectors in radians",
            subtype="ANGLE",
        )
        a_b = tree.outputs.vector(
            "A×B",
            description="Axis around which the angle rotates (cross product of A and B)",
        )

        vector_math = a.normalize()
        vector_math_1 = b.normalize()
        vector_math.dot(vector_math_1).acos() >> angle
        vector_math.cross(vector_math_1) >> a_b


ASSET = VectorAngle

ASSET_METADATA = {
    "description": "The angle between two vectors, in radians",
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
