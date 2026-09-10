# Node-group asset 'Transform Local' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    MatrixSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputMatrix, InputVector


class TransformLocal(AssetGeometryGroup):
    """
    Apply the transform after first moving to world origin, then returning to original position

    Parameters
    ----------
    origin : InputVector
        Vector that defines the local space origin, defaults to `Position`
    transform : InputMatrix
        Transform to apply in local space

    Inputs
    ------
    i.origin : VectorSocket
        Vector that defines the local space origin, defaults to `Position`
    i.transform : MatrixSocket
        Transform to apply in local space

    Outputs
    -------
    o.transform : MatrixSocket
        The final transform
    """

    _name = "Transform Local"
    _asset_name = "Transform Local"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {
        "description": "Apply the transform after first moving to world origin, then returning to original position",
        "node_tool_idname": "geometry.transform_local",
    }

    class _Inputs(SocketAccessor):
        origin: VectorSocket
        """Vector that defines the local space origin, defaults to `Position`"""
        transform: MatrixSocket
        """Transform to apply in local space"""

    class _Outputs(SocketAccessor):
        transform: MatrixSocket
        """The final transform"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        origin: InputVector = None,
        transform: InputMatrix = None,
    ):
        super().__init__(**{"Origin": origin, "Transform": transform})

    def _build_group(self, tree):
        origin = tree.inputs.vector(
            "Origin",
            (0.0, 0.0, 0.0),
            description="Vector that defines the local space origin, defaults to `Position`",
            subtype="TRANSLATION",
            default_input="POSITION",
        )
        transform = tree.inputs.matrix(
            "Transform", description="Transform to apply in local space"
        )
        transform_1 = tree.outputs.matrix(
            "Transform", description="The final transform"
        )

        combine_transform = g.CombineTransform(translation=origin)
        multiply_matrices = g.MultiplyMatrices(
            matrix=g.MultiplyMatrices(matrix=combine_transform, matrix_001=transform),
            matrix_001=combine_transform.o.transform.invert(),
        )

        multiply_matrices >> transform_1


ASSET = TransformLocal

ASSET_METADATA = {
    "description": "Apply the transform after first moving to world origin, then returning to original position",
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
