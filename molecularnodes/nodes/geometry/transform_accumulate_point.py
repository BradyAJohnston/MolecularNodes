# Node-group asset 'Transform Accumulate Point' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    IntegerSocket,
    MatrixSocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import (
    InputBoolean,
    InputInteger,
    InputMatrix,
    InputMenu,
    InputVector,
)
from .transform_accumulate import TransformAccumulate


class TransformAccumulatePoint(AssetGeometryGroup):
    """
    Accumulate transforms on the point domain

    Parameters
    ----------
    domain : InputMenu | Literal["Point", "Edge", "Face", "Face Corner", "Spline", "Instance"]
        Domain on which to accumulate the transforms
    accumulate : InputBoolean
        Include the transform in the final accumlation
    position : InputVector
        Point to transform, defaults to `Position`
    transform : InputMatrix
        Transform field to accumulate
    group_id : InputInteger
        Transform field is accumulated individually for each `Group ID`

    Inputs
    ------
    i.domain : MenuSocket
        Domain on which to accumulate the transforms
    i.accumulate : BooleanSocket
        Include the transform in the final accumlation
    i.position : VectorSocket
        Point to transform, defaults to `Position`
    i.transform : MatrixSocket
        Transform field to accumulate
    i.group_id : IntegerSocket
        Transform field is accumulated individually for each `Group ID`

    Outputs
    -------
    o.vector : VectorSocket
        Transformed vector
    """

    _name = "Transform Accumulate Point"
    _asset_name = "Transform Accumulate Point"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {
        "description": "Accumulate transforms on the point domain",
        "node_tool_idname": "geometry.transform_accumulate_point",
    }

    class _Inputs(SocketAccessor):
        domain: MenuSocket
        """Domain on which to accumulate the transforms"""
        accumulate: BooleanSocket
        """Include the transform in the final accumlation"""
        position: VectorSocket
        """Point to transform, defaults to `Position`"""
        transform: MatrixSocket
        """Transform field to accumulate"""
        group_id: IntegerSocket
        """Transform field is accumulated individually for each `Group ID`"""

    class _Outputs(SocketAccessor):
        vector: VectorSocket
        """Transformed vector"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        domain: InputMenu
        | Literal[
            "Point", "Edge", "Face", "Face Corner", "Spline", "Instance"
        ] = "Point",
        accumulate: InputBoolean = True,
        position: InputVector = None,
        transform: InputMatrix = None,
        group_id: InputInteger = 0,
    ):
        super().__init__(
            **{
                "Domain": domain,
                "Accumulate": accumulate,
                "Position": position,
                "Transform": transform,
                "Group ID": group_id,
            }
        )

    def _build_group(self, tree):
        domain = tree.inputs.menu(
            "Domain",
            description="Domain on which to accumulate the transforms",
            optional_label=True,
        )
        accumulate = tree.inputs.boolean(
            "Accumulate",
            True,
            description="Include the transform in the final accumlation",
            hide_value=True,
        )
        position = tree.inputs.vector(
            "Position",
            (0.0, 0.0, 0.0),
            description="Point to transform, defaults to `Position`",
            hide_value=True,
            subtype="XYZ",
            default_input="POSITION",
        )
        transform = tree.inputs.matrix(
            "Transform", description="Transform field to accumulate"
        )
        group_id = tree.inputs.integer(
            "Group ID",
            0,
            description="Transform field is accumulated individually for each `Group ID`",
            hide_value=True,
        )
        vector = tree.outputs.vector(
            "Vector", description="Transformed vector", subtype="XYZ"
        )

        group = TransformAccumulate(
            domain=domain, accumulate=accumulate, transform=transform, group_id=group_id
        )
        position.transform(group) >> vector

        domain.default_value = "Point"


ASSET = TransformAccumulatePoint

ASSET_METADATA = {
    "description": "Accumulate transforms on the point domain",
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
