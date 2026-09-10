# Node-group asset 'Transform Local Axis' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    MatrixSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputFloat, InputVector


class TransformLocalAxis(AssetGeometryGroup):
    """
    Create a transform around an axis in local space defined by the `Origin` point

    Parameters
    ----------
    origin : InputVector
        The vector defining the local space. Defaults to `Position`
    axis : InputVector
        The axis to rotate around
    angle : InputFloat
        Amount to rotate around the axis

    Inputs
    ------
    i.origin : VectorSocket
        The vector defining the local space. Defaults to `Position`
    i.axis : VectorSocket
        The axis to rotate around
    i.angle : FloatSocket
        Amount to rotate around the axis

    Outputs
    -------
    o.transform : MatrixSocket
        The Transform around the axis
    """

    _name = "Transform Local Axis"
    _asset_name = "Transform Local Axis"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {
        "description": "Create a transform around an axis in local space defined by the `Origin` point",
        "node_tool_idname": "geometry.transform_local_axis",
    }

    class _Inputs(SocketAccessor):
        origin: VectorSocket
        """The vector defining the local space. Defaults to `Position`"""
        axis: VectorSocket
        """The axis to rotate around"""
        angle: FloatSocket
        """Amount to rotate around the axis"""

    class _Outputs(SocketAccessor):
        transform: MatrixSocket
        """The Transform around the axis"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        origin: InputVector = None,
        axis: InputVector = None,
        angle: InputFloat = 0.0,
    ):
        super().__init__(**{"Origin": origin, "Axis": axis, "Angle": angle})

    def _build_group(self, tree):
        origin = tree.inputs.vector(
            "Origin",
            (0.0, 0.0, 0.0),
            description="The vector defining the local space. Defaults to `Position`",
            hide_value=True,
            subtype="TRANSLATION",
            default_input="POSITION",
        )
        axis = tree.inputs.vector(
            "Axis", (0.0, 0.0, 1.0), description="The axis to rotate around"
        )
        angle = tree.inputs.float(
            "Angle",
            0.0,
            description="Amount to rotate around the axis",
            subtype="ANGLE",
        )
        transform = tree.outputs.matrix(
            "Transform", description="The Transform around the axis"
        )

        combine_transform = g.CombineTransform(
            translation=origin, rotation=g.AxisAngleToRotation(axis=axis, angle=angle)
        )
        multiply_matrices = g.MultiplyMatrices(
            matrix=combine_transform,
            matrix_001=g.CombineTransform(translation=origin * -1.0),
        )

        multiply_matrices >> transform


ASSET = TransformLocalAxis

ASSET_METADATA = {
    "description": "Create a transform around an axis in local space defined by the `Origin` point",
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
