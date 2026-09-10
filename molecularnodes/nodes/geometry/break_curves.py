# Node-group asset 'Break Curves' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    GeometrySocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputFloat, InputGeometry
from .edge_length import EdgeLength


class BreakCurves(AssetGeometryGroup):
    """
    Break Curves

    Parameters
    ----------
    curves : InputGeometry
        Geometry to delete elements from
    threshold : InputFloat
        Threshold

    Inputs
    ------
    i.curves : GeometrySocket
        Geometry to delete elements from
    i.threshold : FloatSocket
        Threshold

    Outputs
    -------
    o.curves : GeometrySocket
        Curves
    o.original_index : IntegerSocket
        Index of Spline before segments removed
    """

    _name = "Break Curves"
    _asset_name = "Break Curves"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        curves: GeometrySocket
        """Geometry to delete elements from"""
        threshold: FloatSocket
        """Threshold"""

    class _Outputs(SocketAccessor):
        curves: GeometrySocket
        """Curves"""
        original_index: IntegerSocket
        """Index of Spline before segments removed"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        curves: InputGeometry = None,
        threshold: InputFloat = 0.5,
    ):
        super().__init__(**{"Curves": curves, "Threshold": threshold})

    def _build_group(self, tree):
        curves = tree.inputs.geometry(
            "Curves", description="Geometry to delete elements from"
        )
        threshold = tree.inputs.float("Threshold", 0.5)
        curves_1 = tree.outputs.geometry("Curves")
        original_index = tree.outputs.integer(
            "Original Index", description="Index of Spline before segments removed"
        )

        capture = g.CaptureAttribute.curve(geometry=curves)
        index = capture.items.integer("Index", g.Index())
        (
            capture.o.geometry
            >> g.CurveToMesh()
            >> g.DeleteGeometry.edge(selection=EdgeLength() > threshold)
            >> g.MeshToCurve()
            >> curves_1
        )

        index.output >> original_index


ASSET = BreakCurves

ASSET_METADATA = {
    "catalog_id": "9c167a5c-d0a6-457d-9e9c-90f3edd29e10",
}
