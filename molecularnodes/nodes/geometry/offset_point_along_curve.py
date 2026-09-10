# Node-group asset 'Offset Point Along Curve' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputFloat, InputInteger
from ._shared.index_mixed import IndexMixed
from .between_float import BetweenFloat
from .index_mix_float import IndexMixFloat


class OffsetPointAlongCurve(AssetGeometryGroup):
    """
    Offset along the current point's curve, by a number of points. 1 offsets by a single point on the curve (regardless of how far away they are). 1.5 offsets to half way between `Point Index` + 1 and `Point Index` + 2, returning the `Length` and the `Factor` for this point on the curve

    Parameters
    ----------
    point_index : InputInteger
        The field to evaluate at the given `Index` + `Offset` on the point domain
    offset : InputFloat
        The offset to apply to the `Index` before evaluating the input field

    Inputs
    ------
    i.point_index : IntegerSocket
        The field to evaluate at the given `Index` + `Offset` on the point domain
    i.offset : FloatSocket
        The offset to apply to the `Index` before evaluating the input field

    Outputs
    -------
    o.is_off_spline : BooleanSocket
        The field evaluated at the offset `Index` value
    o.factor : FloatSocket
        The field evaluated at the offset `Index` value
    o.length : FloatSocket
        The field evaluated at the offset `Index` value
    o.index_a : IntegerSocket
        The field evaluated at the offset `Index` value
    o.index_b : IntegerSocket
        The field evaluated at the offset `Index` value
    """

    _name = "Offset Point Along Curve"
    _asset_name = "Offset Point Along Curve"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {
        "description": "Offset along the current point's curve, by a number of points. 1 offsets by a single point on the curve (regardless of how far away they are). 1.5 offsets to half way between `Point Index` + 1 and `Point Index` + 2, returning the `Length` and the `Factor` for this point on the curve",
        "node_tool_idname": "geometry.offset_point_along_curve",
    }

    class _Inputs(SocketAccessor):
        point_index: IntegerSocket
        """The field to evaluate at the given `Index` + `Offset` on the point domain"""
        offset: FloatSocket
        """The offset to apply to the `Index` before evaluating the input field"""

    class _Outputs(SocketAccessor):
        is_off_spline: BooleanSocket
        """The field evaluated at the offset `Index` value"""
        factor: FloatSocket
        """The field evaluated at the offset `Index` value"""
        length: FloatSocket
        """The field evaluated at the offset `Index` value"""
        index_a: IntegerSocket
        """The field evaluated at the offset `Index` value"""
        index_b: IntegerSocket
        """The field evaluated at the offset `Index` value"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        point_index: InputInteger = 0,
        offset: InputFloat = 0.0,
    ):
        super().__init__(**{"Point Index": point_index, "Offset": offset})

    def _build_group(self, tree):
        point_index = tree.inputs.integer(
            "Point Index",
            0,
            description="The field to evaluate at the given `Index` + `Offset` on the point domain",
            hide_value=True,
            default_input="INDEX",
        )
        offset = tree.inputs.float(
            "Offset",
            0.0,
            description="The offset to apply to the `Index` before evaluating the input field",
            min_value=-10_000.0,
            max_value=10_000.0,
        )
        is_off_spline = tree.outputs.boolean(
            "Is Off Spline",
            description="The field evaluated at the offset `Index` value",
        )
        factor = tree.outputs.float(
            "Factor", description="The field evaluated at the offset `Index` value"
        )
        length = tree.outputs.float(
            "Length", description="The field evaluated at the offset `Index` value"
        )
        index_a = tree.outputs.integer(
            "Index A", description="The field evaluated at the offset `Index` value"
        )
        index_b = tree.outputs.integer(
            "Index B", description="The field evaluated at the offset `Index` value"
        )

        with g.Frame("Clamp the offset to not exceed the points in the spline"):
            curve_of_point = g.CurveOfPoint(point_index=point_index)
            integer_math = -curve_of_point.o.index_in_curve
            integer_math_1 = (
                g.PointsOfCurve(curve_index=curve_of_point.o.curve_index).o.total
                - 1
                - curve_of_point.o.index_in_curve
            )
            clamp = offset.clamp(integer_math, integer_math_1)
        group = IndexMixed(index=point_index, offset=clamp)
        (
            BetweenFloat(value=offset, lower=integer_math, upper=integer_math_1)
            >> is_off_spline
        )
        spline_parameter = g.SplineParameter()
        IndexMixFloat(value=spline_parameter.o.factor, index=group.o.mixed) >> factor
        IndexMixFloat(value=spline_parameter.o.length, index=group.o.mixed) >> length

        group.o.floor >> index_a
        group.o.ceiling >> index_b


ASSET = OffsetPointAlongCurve

ASSET_METADATA = {
    "description": "Offset along the current point's curve, by a number of points. 1 offsets by a single point on the curve (regardless of how far away they are). 1.5 offsets to half way between `Point Index` + 1 and `Point Index` + 2, returning the `Length` and the `Factor` for this point on the curve",
    "catalog_id": "9c167a5c-d0a6-457d-9e9c-90f3edd29e10",
}
