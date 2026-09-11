# Node-group asset "Curve Endpoint Values" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputInteger


class CurveEndpointValues(AssetGeometryGroup):
    """
    Output a different integer value for the endpoints of a curve and the middle of a curve

    Parameters
    ----------
    start_size : InputInteger
        The size of the starting points of the curve
    start_value : InputInteger
        The value to output for the starting points of the curve
    other_value : InputInteger
        The value for the points which aren't part of the start or end points
    end_size : InputInteger
        The size of the end points
    end_value : InputInteger
        The value that is returned for the end points of the curve

    Inputs
    ------
    i.start_size : IntegerSocket
        The size of the starting points of the curve
    i.start_value : IntegerSocket
        The value to output for the starting points of the curve
    i.other_value : IntegerSocket
        The value for the points which aren't part of the start or end points
    i.end_size : IntegerSocket
        The size of the end points
    i.end_value : IntegerSocket
        The value that is returned for the end points of the curve

    Outputs
    -------
    o.value : IntegerSocket
        The value for the point, determined by the inputs if it is in the start, end or 'other' region
    o.start_selection : BooleanSocket
        The `n` points at the start of a spline determined by `Start Size`
    o.end_selection : BooleanSocket
        The `n` points at the end of a spline determined by `End Size`
    """

    _name = "Curve Endpoint Values"
    _asset_name = "Curve Endpoint Values"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {
        "description": "Output a different integer value for the endpoints of a curve and the middle of a curve",
        "node_tool_idname": "geometry.curve_endpoint_values",
    }

    class _Inputs(SocketAccessor):
        start_size: IntegerSocket
        """The size of the starting points of the curve"""
        start_value: IntegerSocket
        """The value to output for the starting points of the curve"""
        other_value: IntegerSocket
        """The value for the points which aren't part of the start or end points"""
        end_size: IntegerSocket
        """The size of the end points"""
        end_value: IntegerSocket
        """The value that is returned for the end points of the curve"""

    class _Outputs(SocketAccessor):
        value: IntegerSocket
        """The value for the point, determined by the inputs if it is in the start, end or 'other' region"""
        start_selection: BooleanSocket
        """The `n` points at the start of a spline determined by `Start Size`"""
        end_selection: BooleanSocket
        """The `n` points at the end of a spline determined by `End Size`"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        start_size: InputInteger = 1,
        start_value: InputInteger = 1,
        other_value: InputInteger = 0,
        end_size: InputInteger = 1,
        end_value: InputInteger = -1,
    ):
        super().__init__(
            **{
                "Start Size": start_size,
                "Start Value": start_value,
                "Other Value": other_value,
                "End Size": end_size,
                "End Value": end_value,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        start_size = tree.inputs.integer(
            "Start Size",
            1,
            description="The size of the starting points of the curve",
            min_value=0,
        )
        start_value = tree.inputs.integer(
            "Start Value",
            1,
            description="The value to output for the starting points of the curve",
            min_value=-2147483647,
        )
        other_value = tree.inputs.integer(
            "Other Value",
            0,
            description="The value for the points which aren't part of the start or end points",
            min_value=-2147483647,
        )
        end_size = tree.inputs.integer(
            "End Size", 1, description="The size of the end points", min_value=0
        )
        end_value = tree.inputs.integer(
            "End Value",
            -1,
            description="The value that is returned for the end points of the curve",
        )
        value = tree.outputs.integer(
            "Value",
            description="The value for the point, determined by the inputs if it is in the start, end or 'other' region",
        )
        start_selection = tree.outputs.boolean(
            "Start Selection",
            description="The `n` points at the start of a spline determined by `Start Size`",
        )
        end_selection = tree.outputs.boolean(
            "End Selection",
            description="The `n` points at the end of a spline determined by `End Size`",
        )

        endpoint_selection = g.EndpointSelection(start_size=start_size, end_size=0)
        endpoint_selection_1 = g.EndpointSelection(end_size=end_size, start_size=0)
        (
            endpoint_selection_1.o.selection.switch.integer(
                endpoint_selection.o.selection.switch.integer(other_value, start_value),
                end_value,
            )
            >> value
        )

        endpoint_selection >> start_selection
        endpoint_selection_1 >> end_selection


ASSET = CurveEndpointValues

ASSET_METADATA = {
    "description": "Output a different integer value for the endpoints of a curve and the middle of a curve",
    "catalog_id": "9c167a5c-d0a6-457d-9e9c-90f3edd29e10",
}
