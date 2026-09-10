# Node-group asset 'Normalize Field' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputFloat, InputInteger


class NormalizeField(AssetGeometryGroup):
    """
    Normalize Field

    Parameters
    ----------
    field : InputFloat
        Field
    group_id : InputInteger
        An index used to group values together for multiple separate operations

    Inputs
    ------
    i.field : FloatSocket
        Field
    i.group_id : IntegerSocket
        An index used to group values together for multiple separate operations

    Outputs
    -------
    o.result : FloatSocket
        Result
    """

    _name = "Normalize Field"
    _asset_name = "Normalize Field"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"

    class _Inputs(SocketAccessor):
        field: FloatSocket
        """Field"""
        group_id: IntegerSocket
        """An index used to group values together for multiple separate operations"""

    class _Outputs(SocketAccessor):
        result: FloatSocket
        """Result"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        field: InputFloat = 0.0,
        group_id: InputInteger = 0,
    ):
        super().__init__(**{"Field": field, "Group ID": group_id})

    def _build_group(self, tree):
        field = tree.inputs.float("Field", 0.0)
        group_id = tree.inputs.integer(
            "Group ID",
            0,
            description="An index used to group values together for multiple separate operations",
            hide_value=True,
        )
        result = tree.outputs.float("Result")

        field_min_max = g.FieldMinAndMax.point.float(field, group_id)
        field.map_range(field_min_max.o.min, field_min_max.o.max) >> result


ASSET = NormalizeField

ASSET_METADATA = {
    "description": "Remap a field to between 0 and 1 from it's min and max values.",
    "catalog_id": "7ccb8802-a69f-483e-bf6e-4a47aaa9e940",
}
