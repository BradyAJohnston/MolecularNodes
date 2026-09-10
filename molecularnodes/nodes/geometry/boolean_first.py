# Node-group asset "Boolean First" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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


class BooleanFirst(AssetGeometryGroup):
    """
    Only the first `True` in each `Group ID` remains `True`, all others become `False`

    Parameters
    ----------
    boolean : InputBoolean
        The `Boolean` field to test
    group_id : InputInteger
        Each `Group ID` to find the first `True` element for

    Inputs
    ------
    i.boolean : BooleanSocket
        The `Boolean` field to test
    i.group_id : IntegerSocket
        Each `Group ID` to find the first `True` element for

    Outputs
    -------
    o.is_first : BooleanSocket
        `True` for the first true element in each `Group ID`
    """

    _name = "Boolean First"
    _asset_name = "Boolean First"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {
        "description": "Only the first `True` in each `Group ID` remains `True`, all others become `False`",
        "node_tool_idname": "geometry.boolean_first",
    }

    class _Inputs(SocketAccessor):
        boolean: BooleanSocket
        """The `Boolean` field to test"""
        group_id: IntegerSocket
        """Each `Group ID` to find the first `True` element for"""

    class _Outputs(SocketAccessor):
        is_first: BooleanSocket
        """`True` for the first true element in each `Group ID`"""

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

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        boolean = tree.inputs.boolean(
            "Boolean", False, description="The `Boolean` field to test", hide_value=True
        )
        group_id = tree.inputs.integer(
            "Group ID",
            0,
            description="Each `Group ID` to find the first `True` element for",
            hide_value=True,
        )
        is_first = tree.outputs.boolean(
            "Is First",
            description="`True` for the first true element in each `Group ID`",
        )

        (
            g.Compare.integer.equal(
                g.AccumulateField.point.integer(
                    boolean, group_id
                ).o.leading.point.leading(group_id),
                1,
            )
            >> is_first
        )


ASSET = BooleanFirst

ASSET_METADATA = {
    "description": "Only the first `True` in each `Group ID` remains `True`, all others become `False`",
    "catalog_id": "7ccb8802-a69f-483e-bf6e-4a47aaa9e940",
}
