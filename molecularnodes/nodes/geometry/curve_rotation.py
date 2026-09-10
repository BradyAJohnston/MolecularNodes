# Node-group asset 'Curve Rotation' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    PackageLibrary,
    RotationSocket,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputVector


class CurveRotation(AssetGeometryGroup):
    """
    Computes the rotation of the point on the curve, by default using the `Normal` attribute and the `Tangent` axis

    Parameters
    ----------
    normal : InputVector
        The default direction to use as the Secondary Axis / X axis when computing the rotation

    Inputs
    ------
    i.normal : VectorSocket
        The default direction to use as the Secondary Axis / X axis when computing the rotation

    Outputs
    -------
    o.rotation : RotationSocket
        The computed rotation for the point on the curve
    """

    _name = "Curve Rotation"
    _asset_name = "Curve Rotation"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {
        "description": "Computes the rotation of the point on the curve, by default using the `Normal` attribute and the `Tangent` axis",
        "node_tool_idname": "geometry.curve_rotation",
    }

    class _Inputs(SocketAccessor):
        normal: VectorSocket
        """The default direction to use as the Secondary Axis / X axis when computing the rotation"""

    class _Outputs(SocketAccessor):
        rotation: RotationSocket
        """The computed rotation for the point on the curve"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        normal: InputVector = None,
    ):
        super().__init__(**{"Normal": normal})

    def _build_group(self, tree):
        normal = tree.inputs.vector(
            "Normal",
            (0.0, 0.0, 1.0),
            description="The default direction to use as the Secondary Axis / X axis when computing the rotation",
            default_input="NORMAL",
        )
        rotation = tree.outputs.rotation(
            "Rotation", description="The computed rotation for the point on the curve"
        )

        axes_to_rotation = g.AxesToRotation(
            primary_axis=g.CurveTangent(), secondary_axis=normal
        )

        axes_to_rotation >> rotation


ASSET = CurveRotation

ASSET_METADATA = {
    "description": "Computes the rotation of the point on the curve, by default using the `Normal` attribute and the `Tangent` axis",
    "catalog_id": "9c167a5c-d0a6-457d-9e9c-90f3edd29e10",
}
