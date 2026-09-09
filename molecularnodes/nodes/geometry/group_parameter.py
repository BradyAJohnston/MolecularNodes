# Node-group asset 'Group Parameter' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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


class GroupParameter(AssetGeometryGroup):
    """
    Group Parameter

    Parameters
    ----------
    group_id : InputInteger
        The identifier specifying groupings of the points

    Inputs
    ------
    i.group_id : IntegerSocket
        The identifier specifying groupings of the points

    Outputs
    -------
    o.is_first : BooleanSocket
        If the point is the first point in the `Group ID`
    o.is_last : BooleanSocket
        If the point is the last item in the `Group ID`
    o.group_size : IntegerSocket
        Group Size
    o.relative_index : IntegerSocket
        The relative index of the point within the `Group ID`. Starts at `0` for the first point counting up to `Group Size - 1`
    """

    _name = "Group Parameter"
    _asset_name = "Group Parameter"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"

    class _Inputs(SocketAccessor):
        group_id: IntegerSocket
        """The identifier specifying groupings of the points"""

    class _Outputs(SocketAccessor):
        is_first: BooleanSocket
        """If the point is the first point in the `Group ID`"""
        is_last: BooleanSocket
        """If the point is the last item in the `Group ID`"""
        group_size: IntegerSocket
        """Group Size"""
        relative_index: IntegerSocket
        """The relative index of the point within the `Group ID`. Starts at `0` for the first point counting up to `Group Size - 1`"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        group_id: InputInteger = 0,
    ):
        super().__init__(**{"Group ID": group_id})

    def _build_group(self, tree):
        group_id = tree.inputs.integer(
            "Group ID",
            0,
            description="The identifier specifying groupings of the points",
            hide_value=True,
        )
        is_first = tree.outputs.boolean(
            "Is First", description="If the point is the first point in the `Group ID`"
        )
        is_last = tree.outputs.boolean(
            "Is Last", description="If the point is the last item in the `Group ID`"
        )
        group_size = tree.outputs.integer("Group Size")
        relative_index = tree.outputs.integer(
            "Relative Index",
            description="The relative index  of the point within the `Group ID`. Starts at `0` for the first point counting up to `Group Size - 1`",
        )

        accumulate_field = g.AccumulateField.point.integer(group_index=group_id)
        ~accumulate_field.o.trailing >> is_first
        (
            g.Compare.integer.equal(
                accumulate_field.o.leading, accumulate_field.o.total
            )
            >> is_last
        )

        accumulate_field.o.total >> group_size
        accumulate_field.o.trailing >> relative_index


ASSET = GroupParameter

ASSET_METADATA = {
    "catalog_id": "7ccb8802-a69f-483e-bf6e-4a47aaa9e940",
}
