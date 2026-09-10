# Node-group asset 'Curve Transform' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    MatrixSocket,
    PackageLibrary,
    SocketAccessor,
)
from .curve_rotation import CurveRotation


class CurveTransform(AssetGeometryGroup):
    """
    Calculates the transformation matrix for the point on the curve. Position is taken from the `Position`, `Rotation` is calculated from the `Normal` and `Tangent` values, and the `Radius` drives the scale

    Outputs
    -------
    o.transform : MatrixSocket
        The combine 4X4 transformation matrix for the point of the curve
    """

    _name = "Curve Transform"
    _asset_name = "Curve Transform"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {
        "description": "Calculates the transformation matrix for the point on the curve. Position is taken from the `Position`, `Rotation` is calculated from the `Normal` and `Tangent` values, and the `Radius` drives the scale",
        "node_tool_idname": "geometry.curve_transform",
    }

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        transform: MatrixSocket
        """The combine 4X4 transformation matrix for the point of the curve"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self):
        super().__init__()

    def _build_group(self, tree):
        transform = tree.outputs.matrix(
            "Transform",
            description="The combine 4X4 transformation matrix for the point of the curve",
        )

        combine_transform = g.CombineTransform(
            translation=g.Position(), rotation=CurveRotation(), scale=g.Radius()
        )

        combine_transform >> transform


ASSET = CurveTransform

ASSET_METADATA = {
    "description": "Calculates the transformation matrix for the point on the curve. Position is taken from the `Position`, `Rotation` is calculated from the `Normal` and `Tangent` values, and the `Radius` drives the scale",
    "catalog_id": "9c167a5c-d0a6-457d-9e9c-90f3edd29e10",
}
