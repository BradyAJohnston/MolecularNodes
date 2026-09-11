# Node-group asset "Field Remap" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputFloat, InputInteger


class FieldRemap(AssetGeometryGroup):
    """
    Maps the range of values of the attribute on from the target atoms, to the range from min to max

    Parameters
    ----------
    field : InputFloat
        The field to summarise and remap
    value_min : InputFloat
        The new minimum value of the field
    value_max : InputFloat
        The new maximum value of the field
    group_id : InputInteger
        An index used to group values together for multiple separate operations

    Inputs
    ------
    i.field : FloatSocket
        The field to summarise and remap
    i.value_min : FloatSocket
        The new minimum value of the field
    i.value_max : FloatSocket
        The new maximum value of the field
    i.group_id : IntegerSocket
        An index used to group values together for multiple separate operations

    Outputs
    -------
    o.value : FloatSocket
        The remapped values between the new min and max
    """

    _name = "Field Remap"
    _asset_name = "Field Remap"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {
        "description": "Maps the range of values of the attribute on from the target atoms, to the range from min to max",
        "node_tool_idname": "geometry.field_remap",
    }

    class _Inputs(SocketAccessor):
        field: FloatSocket
        """The field to summarise and remap"""
        value_min: FloatSocket
        """The new minimum value of the field"""
        value_max: FloatSocket
        """The new maximum value of the field"""
        group_id: IntegerSocket
        """An index used to group values together for multiple separate operations"""

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
        field: InputFloat = 0.0,
        value_min: InputFloat = 0.0,
        value_max: InputFloat = 1.0,
        group_id: InputInteger = 0,
    ):
        super().__init__(
            **{
                "Field": field,
                "Value Min": value_min,
                "Value Max": value_max,
                "Group ID": group_id,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
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
        group_id = tree.inputs.integer(
            "Group ID",
            0,
            description="An index used to group values together for multiple separate operations",
            hide_value=True,
        )
        value = tree.outputs.float(
            "Value", description="The remapped values between the new min and max"
        )

        field_min_max = g.FieldMinAndMax.point.float(field, group_id)
        (
            field.map_range(
                field_min_max.o.min, field_min_max.o.max, value_min, value_max
            )
            >> value
        )


ASSET = FieldRemap

ASSET_METADATA = {
    "description": "Maps the range of values of the attribute on from the target atoms, to the range from min to max",
    "catalog_id": "7ccb8802-a69f-483e-bf6e-4a47aaa9e940",
}
