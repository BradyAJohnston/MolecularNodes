# Node-group asset "Select Res ID" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from .boolean_andor import BooleanAndOr


class SelectResID(AssetGeometryGroup):
    """
    Select Res ID

    Parameters
    ----------
    and_ : InputBoolean
        The resulting selection must overlap with this input selection
    or_ : InputBoolean
        The resulting selection can be calculated from this node or be from this input selection
    res_id : InputInteger
        A single `res_id` selection

    Inputs
    ------
    i.and_ : BooleanSocket
        The resulting selection must overlap with this input selection
    i.or_ : BooleanSocket
        The resulting selection can be calculated from this node or be from this input selection
    i.res_id : IntegerSocket
        A single `res_id` selection

    Outputs
    -------
    o.selection : BooleanSocket
        The calculated selection
    o.inverted : BooleanSocket
        The inverse of the calculated selection
    """

    _name = "Select Res ID"
    _asset_name = "Select Res ID"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.select_res_id"}

    class _Inputs(SocketAccessor):
        and_: BooleanSocket
        """The resulting selection must overlap with this input selection"""
        or_: BooleanSocket
        """The resulting selection can be calculated from this node or be from this input selection"""
        res_id: IntegerSocket
        """A single `res_id` selection"""

    class _Outputs(SocketAccessor):
        selection: BooleanSocket
        """The calculated selection"""
        inverted: BooleanSocket
        """The inverse of the calculated selection"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        and_: InputBoolean = True,
        or_: InputBoolean = False,
        res_id: InputInteger = 10,
    ):
        super().__init__(**{"And": and_, "Or": or_, "Res ID": res_id})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        and_ = tree.inputs.boolean(
            "And",
            True,
            description="The resulting selection must overlap with this input selection",
            hide_value=True,
        )
        or_ = tree.inputs.boolean(
            "Or",
            False,
            description="The resulting selection can be calculated from this node or be from this input selection",
            hide_value=True,
        )
        res_id = tree.inputs.integer(
            "Res ID", 10, description="A single `res_id` selection", min_value=0
        )
        selection = tree.outputs.boolean(
            "Selection", description="The calculated selection"
        )
        inverted = tree.outputs.boolean(
            "Inverted", description="The inverse of the calculated selection"
        )

        group = BooleanAndOr(
            and_=and_,
            or_=or_,
            boolean=g.Compare.integer.equal(
                res_id, g.NamedAttribute.integer("res_id").o.attribute
            ),
        )

        group >> selection
        group.o.inverted >> inverted


ASSET = SelectResID

ASSET_METADATA = {
    "catalog_id": "bd1f205b-fea5-4700-b2c2-754f3321e969",
}
