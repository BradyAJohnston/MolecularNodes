# Node-group asset "Group Pick First" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from nodebpy.types import InputBoolean, InputInteger
from .boolean_first import BooleanFirst


class GroupPickFirst(AssetGeometryGroup):
    """
    Get the `Index` of the first `True` item for each `Group ID`

    Parameters
    ----------
    pick : InputBoolean
        Index of the first `True` item in this field will be returned
    group_id : InputInteger
        The first `True` item for each `Group ID` will be returned

    Inputs
    ------
    i.pick : BooleanSocket
        Index of the first `True` item in this field will be returned
    i.group_id : IntegerSocket
        The first `True` item for each `Group ID` will be returned

    Outputs
    -------
    o.is_valid : BooleanSocket
        Valid if as least 1 `True` for the `Group ID`
    o.index : IntegerSocket
        Index of first picked item, returns `-1` if nothing picked
    """

    _name = "Group Pick First"
    _asset_name = "Group Pick First"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {
        "description": "Get the `Index` of the first `True` item for each `Group ID`",
        "node_tool_idname": "geometry.group_pick_first",
    }

    class _Inputs(SocketAccessor):
        pick: BooleanSocket
        """Index of the first `True` item in this field will be returned"""
        group_id: IntegerSocket
        """The first `True` item for each `Group ID` will be returned"""

    class _Outputs(SocketAccessor):
        is_valid: BooleanSocket
        """Valid if as least 1 `True` for the `Group ID`"""
        index: IntegerSocket
        """Index of first picked item, returns `-1` if nothing picked"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        pick: InputBoolean = False,
        group_id: InputInteger = 0,
    ):
        super().__init__(**{"Pick": pick, "Group ID": group_id})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        pick = tree.inputs.boolean(
            "Pick",
            False,
            description="Index of the first `True` item in this field will be returned",
            hide_value=True,
        )
        group_id = tree.inputs.integer(
            "Group ID",
            0,
            description="The first `True` item for each `Group ID` will be returned",
            hide_value=True,
        )
        is_valid = tree.outputs.boolean(
            "Is Valid",
            True,
            description="Valid if as least 1 `True` for the `Group ID`",
        )
        index = tree.outputs.integer(
            "Index",
            description="Index of first picked item, returns `-1` if nothing picked",
            min_value=0,
        )

        compare = g.AccumulateField.point.integer(pick, group_id).o.total >= 1
        accumulate_field = (
            BooleanFirst(boolean=pick, group_id=group_id)
            .o.is_first.switch.integer(true=g.Index())
            .point.total(group_id)
        )
        compare.switch.integer(-1, accumulate_field) >> index

        compare >> is_valid


ASSET = GroupPickFirst

ASSET_METADATA = {
    "description": "Get the `Index` of the first `True` item for each `Group ID`",
    "catalog_id": "7ccb8802-a69f-483e-bf6e-4a47aaa9e940",
}
