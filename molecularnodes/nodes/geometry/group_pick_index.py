# Node-group asset 'Group Pick Index' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputInteger
from .group_parameter import GroupParameter
from .group_pick import GroupPick


class GroupPickIndex(AssetGeometryGroup):
    """
    Group Pick Index

    Parameters
    ----------
    relative_index : InputInteger
        Relative Index
    group_id : InputInteger
        Group ID

    Inputs
    ------
    i.relative_index : IntegerSocket
        Relative Index
    i.group_id : IntegerSocket
        Group ID

    Outputs
    -------
    o.is_valid : BooleanSocket
        If the requested `Index in Group` is valid for this `Group ID`
    o.index : IntegerSocket
        The `Index` of the requested point, within the overall geometry
    """

    _name = "Group Pick Index"
    _asset_name = "Group Pick Index"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"

    class _Inputs(SocketAccessor):
        relative_index: IntegerSocket
        """Relative Index"""
        group_id: IntegerSocket
        """Group ID"""

    class _Outputs(SocketAccessor):
        is_valid: BooleanSocket
        """If the requested `Index in Group` is valid for this `Group ID`"""
        index: IntegerSocket
        """The `Index` of the requested point, within the overall geometry"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        relative_index: InputInteger = 0,
        group_id: InputInteger = 0,
    ):
        super().__init__(**{"Relative Index": relative_index, "Group ID": group_id})

    def _build_group(self, tree):
        relative_index = tree.inputs.integer(
            "Relative Index", 0, min_value=0, hide_value=True
        )
        group_id = tree.inputs.integer("Group ID", 0, hide_value=True)
        is_valid = tree.outputs.boolean(
            "Is Valid",
            True,
            description="If the requested `Index in Group` is valid for this `Group ID`",
        )
        index = tree.outputs.integer(
            "Index",
            description="The `Index` of the requested point, within the overall geometry",
            min_value=0,
        )

        compare = g.Compare.integer.equal(
            relative_index, GroupParameter(group_id=group_id).o.relative_index
        )
        group = GroupPick(pick=compare, group_id=group_id)

        group >> is_valid
        group.o.index >> index


ASSET = GroupPickIndex

ASSET_METADATA = {
    "catalog_id": "7ccb8802-a69f-483e-bf6e-4a47aaa9e940",
}
