# Node-group asset "Group Pick" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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


class GroupPick(AssetGeometryGroup):
    """
    Get the item of the `True` in each `Group ID`, but only if there is a single `True` in each group

    Parameters
    ----------
    pick : InputBoolean
        Selection for the point to pick for the `Group ID`. Will only be valid for a single `Pick` for each `Group ID`
    group_id : InputInteger
        Field definining the `Group ID` to pick from for the points

    Inputs
    ------
    i.pick : BooleanSocket
        Selection for the point to pick for the `Group ID`. Will only be valid for a single `Pick` for each `Group ID`
    i.group_id : IntegerSocket
        Field definining the `Group ID` to pick from for the points

    Outputs
    -------
    o.is_valid : BooleanSocket
        Valid for the point only if there is 1 single `True` for the `Pick` field in the point's group in the `Group ID`
    o.index : IntegerSocket
        Index of picked item. Returns -1 if not a valid pick
    """

    _name = "Group Pick"
    _asset_name = "Group Pick"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {
        "description": "Get the item of the `True` in each `Group ID`, but only if there is a single `True` in each group",
        "node_tool_idname": "geometry.group_pick",
    }

    class _Inputs(SocketAccessor):
        pick: BooleanSocket
        """Selection for the point to pick for the `Group ID`. Will only be valid for a single `Pick` for each `Group ID`"""
        group_id: IntegerSocket
        """Field definining the `Group ID` to pick from for the points"""

    class _Outputs(SocketAccessor):
        is_valid: BooleanSocket
        """Valid for the point only if there is 1 single `True` for the `Pick` field in the point's group in the `Group ID`"""
        index: IntegerSocket
        """Index of picked item. Returns -1 if not a valid pick"""

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
            description="Selection for the point to pick for the `Group ID`. Will only be valid for a single `Pick` for each `Group ID`",
            hide_value=True,
        )
        group_id = tree.inputs.integer(
            "Group ID",
            0,
            description="Field definining the `Group ID` to pick from for the points",
            hide_value=True,
        )
        is_valid = tree.outputs.boolean(
            "Is Valid",
            True,
            description="Valid for the point only if there is 1 single `True` for the `Pick` field in the point's group in the `Group ID`",
        )
        index = tree.outputs.integer(
            "Index",
            description="Index of picked item. Returns -1 if not a valid pick",
            min_value=0,
        )

        compare = g.Compare.integer.equal(
            g.AccumulateField.point.integer(pick, group_id).o.total, 1
        )
        (
            compare.o.result.switch.integer(
                -1, pick.switch.integer(true=g.Index()).point.total(group_id)
            )
            >> index
        )

        compare >> is_valid


ASSET = GroupPick

ASSET_METADATA = {
    "description": "Get the item of the `True` in each `Group ID`, but only if there is a single `True` in each group",
    "catalog_id": "7ccb8802-a69f-483e-bf6e-4a47aaa9e940",
}
