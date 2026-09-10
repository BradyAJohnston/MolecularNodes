# Node-group asset "Curve Vectors" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)


class CurveVectors(AssetGeometryGroup):
    """
    Curve Vectors

    Outputs
    -------
    o.normal : VectorSocket
        The normal of the control point. Used for calculating the rotation for the point and when calculating `Curve to Mesh`
    o.tangent : VectorSocket
        The tangent of the point, which is calculated as the direction from the previous point to the next point
    o.bitangent : VectorSocket
        The cross product of the Normal and the Tangent
    """

    _name = "Curve Vectors"
    _asset_name = "Curve Vectors"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.curve_vectors"}

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        normal: VectorSocket
        """The normal of the control point. Used for calculating the rotation for the point and when calculating `Curve to Mesh`"""
        tangent: VectorSocket
        """The tangent of the point, which is calculated as the direction from the previous point to the next point"""
        bitangent: VectorSocket
        """The cross product of the Normal and the Tangent"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self):
        super().__init__()

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        normal = tree.outputs.vector(
            "Normal",
            description="The normal of the control point. Used for calculating the rotation for the point and when calculating `Curve to Mesh`",
        )
        tangent = tree.outputs.vector(
            "Tangent",
            description="The tangent of the point, which is calculated as the direction from the previous point to the next point",
        )
        bitangent = tree.outputs.vector(
            "Bitangent", description="The cross product of the Normal and the Tangent"
        )

        curve_tangent = g.CurveTangent()
        normal_1 = g.Normal(legacy_corner_normals=True)
        curve_tangent.o.tangent.cross(normal_1.o.normal) >> bitangent

        normal_1 >> normal
        curve_tangent >> tangent


ASSET = CurveVectors

ASSET_METADATA = {
    "catalog_id": "9c167a5c-d0a6-457d-9e9c-90f3edd29e10",
}
