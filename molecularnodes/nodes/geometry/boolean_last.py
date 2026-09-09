# Node-group asset 'Boolean Last' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from nodebpy.types import InputBoolean
from .group_info import GroupInfo
from .offset_integer import OffsetInteger


class BooleanLast(AssetGeometryGroup):
    """
    Index of last time the `Boolean` is true for each `Group ID` (not including the current point).

    Parameters
    ----------
    boolean : InputBoolean
        Value to test for True when tracking `Index` locations

    Inputs
    ------
    i.boolean : BooleanSocket
        Value to test for True when tracking `Index` locations

    Outputs
    -------
    o.index_of_last : IntegerSocket
        Accumulating in the point domain, this is the index where the `Boolean` was last `True`. For a point where the `Boolean` is `True`, the index will be of the _previous_ time this was true, but the next point will then reference this point
    """

    _name = "Boolean Last"
    _asset_name = "Boolean Last"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {
        "description": "Index of last time the `Boolean` is true for each `Group ID` (not including the current point). "
    }

    class _Inputs(SocketAccessor):
        boolean: BooleanSocket
        """Value to test for True when tracking `Index` locations"""

    class _Outputs(SocketAccessor):
        index_of_last: IntegerSocket
        """Accumulating in the point domain, this is the index where the `Boolean` was last `True`. For a point where the `Boolean` is `True`, the index will be of the _previous_ time this was true, but the next point will then reference this point"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        boolean: InputBoolean = True,
    ):
        super().__init__(**{"Boolean": boolean})

    def _build_group(self, tree):
        boolean = tree.inputs.boolean(
            "Boolean",
            True,
            description="Value to test for True when tracking `Index` locations",
            hide_value=True,
        )
        index_of_last = tree.outputs.integer(
            "Index of Last",
            description="Accumulating in the point domain, this is the index where the `Boolean` was last `True`. For a point where the `Boolean` is `True`, the index will be of the _previous_ time this was true, but the next point will then reference this point",
        )

        accumulate_field = g.AccumulateField.point.integer(boolean)
        compare = g.Compare.integer.equal(
            g.AccumulateField.point.integer(
                group_index=accumulate_field.o.leading
            ).o.trailing,
            0,
        )
        group = OffsetInteger(
            integer=GroupInfo(group_id=accumulate_field.o.leading).o.index_of_first,
            offset=-1,
        )
        _group_1 = OffsetInteger(
            integer=compare.o.result.switch.integer(true=g.Index()).point.total(
                accumulate_field.o.trailing
            )
        )
        (accumulate_field.o.leading > 0).switch.integer(-1, group) >> index_of_last


ASSET = BooleanLast

ASSET_METADATA = {
    "description": "Index of last time the `Boolean` is true for each `Group ID` (not including the current point). ",
    "catalog_id": "7ccb8802-a69f-483e-bf6e-4a47aaa9e940",
}
