# Node-group asset "Integer Run" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from .offset_integer import OffsetInteger


class IntegerRun(AssetGeometryGroup):
    """
    A unique value for each grouping of a value. Accumulating along the field, the output Group Mask increments by 1 whenever the value or Group ID changes

    Parameters
    ----------
    value : InputInteger
        The field to check for changes in value
    group_id : InputInteger
        Does not restart counting for each Group ID, but does increment by 1 when the Group ID changes

    Inputs
    ------
    i.value : IntegerSocket
        The field to check for changes in value
    i.group_id : IntegerSocket
        Does not restart counting for each Group ID, but does increment by 1 when the Group ID changes

    Outputs
    -------
    o.is_different : BooleanSocket
        The current value is different from the previous
    o.group_id : IntegerSocket
        The new `Group ID`, which increases by 1 whenever the `Value` or `Group ID` values change
    """

    _name = "Integer Run"
    _asset_name = "Integer Run"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {
        "description": "A unique value for each grouping of a value. Accumulating along the field, the output Group Mask increments by 1 whenever the value or Group ID changes",
        "node_tool_idname": "geometry.integer_run",
    }

    class _Inputs(SocketAccessor):
        value: IntegerSocket
        """The field to check for changes in value"""
        group_id: IntegerSocket
        """Does not restart counting for each Group ID, but does increment by 1 when the Group ID changes"""

    class _Outputs(SocketAccessor):
        is_different: BooleanSocket
        """The current value is different from the previous"""
        group_id: IntegerSocket
        """The new `Group ID`, which increases by 1 whenever the `Value` or `Group ID` values change"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        value: InputInteger = 0,
        group_id: InputInteger = 0,
    ):
        super().__init__(**{"Value": value, "Group ID": group_id})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        value = tree.inputs.integer(
            "Value",
            0,
            description="The field to check for changes in value",
            hide_value=True,
        )
        group_id = tree.inputs.integer(
            "Group ID",
            0,
            description="Does not restart counting for each Group ID, but does increment by 1 when the Group ID changes",
            hide_value=True,
        )
        is_different = tree.outputs.boolean(
            "Is Different",
            description="The current value is different from the previous",
        )
        group_id_1 = tree.outputs.integer(
            "Group ID",
            description="The new `Group ID`, which increases by 1 whenever the `Value` or `Group ID` values change",
        )

        boolean_math = g.Compare.integer.not_equal(
            OffsetInteger(integer=value, offset=-1), value
        ).o.result | g.Compare.integer.not_equal(
            group_id, OffsetInteger(integer=group_id, offset=-1)
        )
        accumulate_field = g.AccumulateField.point.integer(boolean_math)

        boolean_math >> is_different
        accumulate_field >> group_id_1


ASSET = IntegerRun

ASSET_METADATA = {
    "description": "A unique value for each grouping of a value. Accumulating along the field, the output Group Mask increments by 1 whenever the value or Group ID changes",
    "catalog_id": "7ccb8802-a69f-483e-bf6e-4a47aaa9e940",
}
