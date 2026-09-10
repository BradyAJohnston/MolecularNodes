# Node-group asset "Offset Curve" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    GeometrySocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputFloat, InputGeometry


class OffsetCurve(AssetGeometryGroup):
    """
    Offset Curve

    Parameters
    ----------
    curve : InputGeometry
        Curve
    points : InputFloat
        Number of points to offset and sample position and normal by. Fractions sample between the points along the curve.

    Inputs
    ------
    i.curve : GeometrySocket
        Curve
    i.points : FloatSocket
        Number of points to offset and sample position and normal by. Fractions sample between the points along the curve.

    Outputs
    -------
    o.curve : GeometrySocket
        Curve
    o.factor : FloatSocket
        The pre-offset factor for the curve's point
    """

    _name = "Offset Curve"
    _asset_name = "Offset Curve"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        curve: GeometrySocket
        """Curve"""
        points: FloatSocket
        """Number of points to offset and sample position and normal by. Fractions sample between the points along the curve."""

    class _Outputs(SocketAccessor):
        curve: GeometrySocket
        """Curve"""
        factor: FloatSocket
        """The pre-offset factor for the curve's point"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        curve: InputGeometry = None,
        points: InputFloat = 0.0,
    ):
        super().__init__(**{"Curve": curve, "Points": points})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        curve = tree.inputs.geometry("Curve")
        points = tree.inputs.float(
            "Points",
            0.0,
            description="Number of points to offset and sample position and normal by. Fractions sample between the points along the curve.",
            min_value=-10_000.0,
            max_value=10_000.0,
        )
        curve_1 = tree.outputs.geometry("Curve")
        factor = tree.outputs.float(
            "Factor", description="The pre-offset factor for the curve's point"
        )

        with g.Frame("Amount of factor for this number of points to offset"):
            math_1 = points / g.SplineLength().o.point_count
        capture = g.CaptureAttribute.point(geometry=curve)
        value = capture.items.float("Value", g.SplineParameter().o.factor + math_1)
        sample_curve = g.SampleCurve(
            curves=capture.o.geometry,
            factor=value.output,
            curve_index=g.CurveOfPoint().o.curve_index,
        )
        capture_1 = g.CaptureAttribute.point(geometry=capture.o.geometry)
        position = capture_1.items.vector("Position", sample_curve.o.position)
        normal = capture_1.items.vector("Normal", sample_curve.o.normal)
        (
            capture_1.o.geometry
            >> g.SetPosition(position=position.output)
            >> g.SetCurveNormal(normal=normal.output, mode="Free")
            >> curve_1
        )

        value.output >> factor


ASSET = OffsetCurve

ASSET_METADATA = {
    "catalog_id": "a1e4128a-131f-4e0e-b54e-81f863aba707",
}
