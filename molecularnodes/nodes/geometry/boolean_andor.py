# Node-group asset 'Boolean AndOr' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean


class BooleanAndOr(AssetGeometryGroup):
    """
    Boolean AndOr

    Parameters
    ----------
    and_ : InputBoolean
        The resulting selection must overlap with this input selection
    or_ : InputBoolean
        The resulting selection can be calculated from this node or be from this input selection
    boolean : InputBoolean
        Boolean

    Inputs
    ------
    i.and_ : BooleanSocket
        The resulting selection must overlap with this input selection
    i.or_ : BooleanSocket
        The resulting selection can be calculated from this node or be from this input selection
    i.boolean : BooleanSocket
        Boolean

    Outputs
    -------
    o.boolean : BooleanSocket
        Boolean
    o.inverted : BooleanSocket
        The inverse of the calculated selection
    """

    _name = "Boolean AndOr"
    _asset_name = "Boolean AndOr"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {"node_tool_idname": "geometry.boolean_andor"}

    class _Inputs(SocketAccessor):
        and_: BooleanSocket
        """The resulting selection must overlap with this input selection"""
        or_: BooleanSocket
        """The resulting selection can be calculated from this node or be from this input selection"""
        boolean: BooleanSocket
        """Boolean"""

    class _Outputs(SocketAccessor):
        boolean: BooleanSocket
        """Boolean"""
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
        boolean: InputBoolean = True,
    ):
        super().__init__(**{"And": and_, "Or": or_, "Boolean": boolean})

    def _build_group(self, tree):
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
        boolean = tree.inputs.boolean("Boolean", True)
        boolean_1 = tree.outputs.boolean("Boolean")
        inverted = tree.outputs.boolean(
            "Inverted", description="The inverse of the calculated selection"
        )

        boolean_math = and_ & boolean | or_
        ~boolean_math >> inverted

        boolean_math >> boolean_1


ASSET = BooleanAndOr

ASSET_METADATA = {
    "catalog_id": "0d42d20f-33f0-464a-b4c4-9fb278826421",
}
