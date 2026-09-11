# Node-group asset "Accumulate Axis Rotation" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    IntegerSocket,
    MatrixSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputBoolean, InputFloat, InputInteger, InputVector
from .boolean_last import BooleanLast
from .transform_accumulate import TransformAccumulate
from .transform_local_axis import TransformLocalAxis


class AccumulateAxisRotation(AssetGeometryGroup):
    """
    Accumulate Axis Rotation

    Parameters
    ----------
    position : InputVector
        Position vector to transform
    selection : InputBoolean
        Selection
    pivot : InputBoolean
        The points at which points where `Accumulate` is true will look to for their axis of transformation.
    angle : InputFloat
        Amount to rotate around the axis
    group_id : InputInteger
        Transform field is accumulated individually for each `Group ID`
    transform_index : InputInteger
        Index at which to evaluate the final transform for the positon. For most cases this will be `Index`, but it might be that some points need to use the accumulated transform from another point instead

    Inputs
    ------
    i.position : VectorSocket
        Position vector to transform
    i.selection : BooleanSocket
        Selection
    i.pivot : BooleanSocket
        The points at which points where `Accumulate` is true will look to for their axis of transformation.
    i.angle : FloatSocket
        Amount to rotate around the axis
    i.group_id : IntegerSocket
        Transform field is accumulated individually for each `Group ID`
    i.transform_index : IntegerSocket
        Index at which to evaluate the final transform for the positon. For most cases this will be `Index`, but it might be that some points need to use the accumulated transform from another point instead

    Outputs
    -------
    o.position : VectorSocket
        Transformed vector
    o.trasnform : MatrixSocket
        The accumlated transform, not yet applied to the `Position` vector
    """

    _name = "Accumulate Axis Rotation"
    _asset_name = "Accumulate Axis Rotation"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"

    class _Inputs(SocketAccessor):
        position: VectorSocket
        """Position vector to transform"""
        selection: BooleanSocket
        """Selection"""
        pivot: BooleanSocket
        """The points at which points where `Accumulate` is true will look to for their axis of transformation."""
        angle: FloatSocket
        """Amount to rotate around the axis"""
        group_id: IntegerSocket
        """Transform field is accumulated individually for each `Group ID`"""
        transform_index: IntegerSocket
        """Index at which to evaluate the final transform for the positon. For most cases this will be `Index`, but it might be that some points need to use the accumulated transform from another point instead"""

    class _Outputs(SocketAccessor):
        position: VectorSocket
        """Transformed vector"""
        trasnform: MatrixSocket
        """The accumlated transform, not yet applied to the `Position` vector"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        position: InputVector = None,
        selection: InputBoolean = True,
        pivot: InputBoolean = False,
        angle: InputFloat = 0.0,
        group_id: InputInteger = 0,
        transform_index: InputInteger = 0,
    ):
        super().__init__(
            **{
                "Position": position,
                "Selection": selection,
                "Pivot": pivot,
                "Angle": angle,
                "Group ID": group_id,
                "Transform Index": transform_index,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        position = tree.inputs.vector(
            "Position",
            (0.0, 0.0, 0.0),
            description="Position vector to transform",
            default_input="POSITION",
        )
        selection = tree.inputs.boolean("Selection", True, hide_value=True)
        pivot = tree.inputs.boolean(
            "Pivot",
            False,
            description="The points at which points where `Accumulate` is true will look to for their axis of transformation.",
            hide_value=True,
        )
        angle = tree.inputs.float(
            "Angle",
            0.0,
            description="Amount to rotate around the axis",
            subtype="ANGLE",
        )
        group_id = tree.inputs.integer(
            "Group ID",
            0,
            description="Transform field is accumulated individually for each `Group ID`",
            hide_value=True,
        )
        transform_index = tree.inputs.integer(
            "Transform Index",
            0,
            description="Index at which to evaluate the final transform for the positon. For most cases this will be `Index`, but it might be that some points need to use the accumulated transform from another point instead",
            min_value=0,
            hide_value=True,
            default_input="INDEX",
        )
        position_1 = tree.outputs.vector(
            "Position", description="Transformed vector", subtype="XYZ"
        )
        trasnform = tree.outputs.matrix(
            "Trasnform",
            description="The accumlated transform, not yet applied to the `Position` vector",
        )

        group = TransformLocalAxis(
            origin=position,
            axis=position.point.at(BooleanLast(boolean=pivot)) - position,
            angle=angle,
        )
        evaluate_at_index = TransformAccumulate(
            accumulate=selection, transform=group, group_id=group_id
        ).o.transform.point.at(transform_index)
        position.transform(evaluate_at_index) >> position_1

        evaluate_at_index >> trasnform


ASSET = AccumulateAxisRotation

ASSET_METADATA = {
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
