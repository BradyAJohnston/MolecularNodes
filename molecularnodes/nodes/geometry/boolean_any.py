# Node-group asset 'Boolean Any' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from nodebpy.types import InputBoolean, InputInteger


class BooleanAny(AssetGeometryGroup):
    """
    Boolean Any

    Parameters
    ----------
    boolean : InputBoolean
        Boolean
    group_id : InputInteger
        Group ID

    Inputs
    ------
    i.boolean : BooleanSocket
        Boolean
    i.group_id : IntegerSocket
        Group ID

    Outputs
    -------
    o.boolean : BooleanSocket
        Boolean
    """

    _name = "Boolean Any"
    _asset_name = "Boolean Any"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"

    class _Inputs(SocketAccessor):
        boolean: BooleanSocket
        """Boolean"""
        group_id: IntegerSocket
        """Group ID"""

    class _Outputs(SocketAccessor):
        boolean: BooleanSocket
        """Boolean"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        boolean: InputBoolean = False,
        group_id: InputInteger = 0,
    ):
        super().__init__(**{"Boolean": boolean, "Group ID": group_id})

    def _build_group(self, tree):
        boolean = tree.inputs.boolean("Boolean", False)
        group_id = tree.inputs.integer("Group ID", 0, hide_value=True)
        boolean_1 = tree.outputs.boolean("Boolean")

        (g.AccumulateField.point.integer(boolean, group_id).o.total > 0) >> boolean_1


ASSET = BooleanAny

ASSET_METADATA = {
    "catalog_id": "7ccb8802-a69f-483e-bf6e-4a47aaa9e940",
}
