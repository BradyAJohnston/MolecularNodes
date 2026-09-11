# Node-group asset "Transform Relative" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    MatrixSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputVector


class TransformRelative(AssetGeometryGroup):
    """
    Transform Relative

    Parameters
    ----------
    a : InputVector
        The final position
    b : InputVector
        The position moving from
    c : InputVector
        The position defining the axis of CB

    Inputs
    ------
    i.a : VectorSocket
        The final position
    i.b : VectorSocket
        The position moving from
    i.c : VectorSocket
        The position defining the axis of CB

    Outputs
    -------
    o.transform : MatrixSocket
        The `Transform` to move from B to A, relative to the axis of C to B
    """

    _name = "Transform Relative"
    _asset_name = "Transform Relative"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {"node_tool_idname": "geometry.transform_relative"}

    class _Inputs(SocketAccessor):
        a: VectorSocket
        """The final position"""
        b: VectorSocket
        """The position moving from"""
        c: VectorSocket
        """The position defining the axis of CB"""

    class _Outputs(SocketAccessor):
        transform: MatrixSocket
        """The `Transform` to move from B to A, relative to the axis of C to B"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        a: InputVector = None,
        b: InputVector = None,
        c: InputVector = None,
    ):
        super().__init__(**{"A": a, "B": b, "C": c})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        a = tree.inputs.vector(
            "A",
            (0.0, 0.0, 0.0),
            description="The final position",
            min_value=-10_000.0,
            max_value=10_000.0,
        )
        b = tree.inputs.vector(
            "B", (0.0, 0.0, 0.0), description="The position moving from"
        )
        c_ = tree.inputs.vector(
            "C",
            (0.0, 0.0, 0.0),
            description="The position defining the axis of CB",
            min_value=-10_000.0,
            max_value=10_000.0,
        )
        transform = tree.outputs.matrix(
            "Transform",
            description="The `Transform` to move from B to A, relative to the axis of C to B",
        )

        vector_math = a - b
        rotate_rotation = g.AlignRotationToVector(vector=vector_math).o.rotation.rotate(
            g.AlignRotationToVector(vector=b - c_).o.rotation.invert()
        )
        combine_transform = g.CombineTransform(
            translation=g.CombineXYZ(z=vector_math.length()), rotation=rotate_rotation
        )

        combine_transform >> transform


ASSET = TransformRelative

ASSET_METADATA = {
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
