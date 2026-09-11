# Node-group asset "Group Pick Vector" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputBoolean, InputInteger, InputVector
from .group_pick import GroupPick


class GroupPickVector(AssetGeometryGroup):
    """
    Group Pick Vector

    Parameters
    ----------
    pick : InputBoolean
        Selection for the point to pick for the `Group ID`. Will only be valid for a single `Pick` for each `Group ID`
    group_id : InputInteger
        Field definining the `Group ID` to pick from for the points
    position : InputVector
        Vector field to pick vlaue for, defaults to `Position`

    Inputs
    ------
    i.pick : BooleanSocket
        Selection for the point to pick for the `Group ID`. Will only be valid for a single `Pick` for each `Group ID`
    i.group_id : IntegerSocket
        Field definining the `Group ID` to pick from for the points
    i.position : VectorSocket
        Vector field to pick vlaue for, defaults to `Position`

    Outputs
    -------
    o.is_valid : BooleanSocket
        Whether the pick for the point's `Group ID` is valid
    o.index : IntegerSocket
        Picked Index for the Group, -1 if not valid
    o.vector : VectorSocket
        Picked vector for the group, `(0, 0, 0)` if not valid
    """

    _name = "Group Pick Vector"
    _asset_name = "Group Pick Vector"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {"node_tool_idname": "geometry.group_pick_vector"}

    class _Inputs(SocketAccessor):
        pick: BooleanSocket
        """Selection for the point to pick for the `Group ID`. Will only be valid for a single `Pick` for each `Group ID`"""
        group_id: IntegerSocket
        """Field definining the `Group ID` to pick from for the points"""
        position: VectorSocket
        """Vector field to pick vlaue for, defaults to `Position`"""

    class _Outputs(SocketAccessor):
        is_valid: BooleanSocket
        """Whether the pick for the point's `Group ID` is valid"""
        index: IntegerSocket
        """Picked Index for the Group, -1 if not valid"""
        vector: VectorSocket
        """Picked vector for the group, `(0, 0, 0)` if not valid"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        pick: InputBoolean = False,
        group_id: InputInteger = 0,
        position: InputVector = None,
    ):
        super().__init__(**{"Pick": pick, "Group ID": group_id, "Position": position})

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
        position = tree.inputs.vector(
            "Position",
            (0.0, 0.0, 0.0),
            description="Vector field to pick vlaue for, defaults to `Position`",
            hide_value=True,
            default_input="POSITION",
        )
        is_valid = tree.outputs.boolean(
            "Is Valid",
            description="Whether the pick for the point's `Group ID` is valid",
        )
        index = tree.outputs.integer(
            "Index", description="Picked Index for the Group, -1 if not valid"
        )
        vector = tree.outputs.vector(
            "Vector",
            description="Picked vector for the group, `(0, 0, 0)` if not valid",
        )

        group = GroupPick(pick=pick, group_id=group_id)
        (
            group.o.is_valid.switch.vector(
                (0.0, 0.0, 0.0), position.point.at(group.o.index)
            )
            >> vector
        )

        group >> is_valid
        group.o.index >> index


ASSET = GroupPickVector

ASSET_METADATA = {
    "catalog_id": "7ccb8802-a69f-483e-bf6e-4a47aaa9e940",
}
