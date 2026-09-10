# Node-group asset "Geometry Field Remap" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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


class GeometryFieldRemap(AssetGeometryGroup):
    """
    Maps the range of values of the attribute on from the target atoms, to the range from min to max

    Parameters
    ----------
    geometry : InputGeometry
        Geometry to get the statistics from
    field : InputFloat
        The field to summarise and remap
    value_min : InputFloat
        The new minimum value of the field
    value_max : InputFloat
        The new maximum value of the field

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry to get the statistics from
    i.field : FloatSocket
        The field to summarise and remap
    i.value_min : FloatSocket
        The new minimum value of the field
    i.value_max : FloatSocket
        The new maximum value of the field

    Outputs
    -------
    o.value : FloatSocket
        The remapped values between the new min and max
    """

    _name = "Geometry Field Remap"
    _asset_name = "Geometry Field Remap"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {
        "description": "Maps the range of values of the attribute on from the target atoms, to the range from min to max",
        "node_tool_idname": "geometry.field_remap",
    }

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry to get the statistics from"""
        field: FloatSocket
        """The field to summarise and remap"""
        value_min: FloatSocket
        """The new minimum value of the field"""
        value_max: FloatSocket
        """The new maximum value of the field"""

    class _Outputs(SocketAccessor):
        value: FloatSocket
        """The remapped values between the new min and max"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        geometry: InputGeometry = None,
        field: InputFloat = 0.0,
        value_min: InputFloat = 0.0,
        value_max: InputFloat = 1.0,
    ):
        super().__init__(
            **{
                "Geometry": geometry,
                "Field": field,
                "Value Min": value_min,
                "Value Max": value_max,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry(
            "Geometry", description="Geometry to get the statistics from"
        )
        field = tree.inputs.float(
            "Field",
            0.0,
            description="The field to summarise and remap",
            hide_value=True,
        )
        value_min = tree.inputs.float(
            "Value Min",
            0.0,
            description="The new minimum value of the field",
            min_value=-10_000.0,
            max_value=10_000.0,
        )
        value_max = tree.inputs.float(
            "Value Max",
            1.0,
            description="The new maximum value of the field",
            min_value=-10_000.0,
            max_value=10_000.0,
        )
        value = tree.outputs.float(
            "Value", description="The remapped values between the new min and max"
        )

        attribute_statistic = g.AttributeStatistic.point.float(
            geometry, attribute=field
        )
        (
            field.map_range(
                attribute_statistic.o.min,
                attribute_statistic.o.max,
                value_min,
                value_max,
            )
            >> value
        )


ASSET = GeometryFieldRemap

ASSET_METADATA = {
    "description": "Maps the range of values of the attribute on from the target atoms, to the range from min to max",
    "catalog_id": "7ccb8802-a69f-483e-bf6e-4a47aaa9e940",
}
