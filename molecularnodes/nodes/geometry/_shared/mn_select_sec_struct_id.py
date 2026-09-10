# Node group '.MN_select_sec_struct_id' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    BooleanSocket,
    CustomGeometryGroup,
    IntegerSocket,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputInteger
from ..boolean_andor import BooleanAndOr


class MN_select_sec_struct_id(CustomGeometryGroup):
    """
    .MN_select_sec_struct_id

    Parameters
    ----------
    and_ : InputBoolean
        The resulting selection must overlap with this input selection
    or_ : InputBoolean
        The resulting selection can be calculated from this node or be from this input selection
    id : InputInteger
        Secondary structure component to select

    Inputs
    ------
    i.and_ : BooleanSocket
        The resulting selection must overlap with this input selection
    i.or_ : BooleanSocket
        The resulting selection can be calculated from this node or be from this input selection
    i.id : IntegerSocket
        Secondary structure component to select

    Outputs
    -------
    o.selection : BooleanSocket
        The calculated selection
    o.inverted : BooleanSocket
        The inverse of the calculated selection
    """

    _name = ".MN_select_sec_struct_id"
    _tree_properties = {"node_tool_idname": "geometry._mn_select_sec_struct_id"}

    class _Inputs(SocketAccessor):
        and_: BooleanSocket
        """The resulting selection must overlap with this input selection"""
        or_: BooleanSocket
        """The resulting selection can be calculated from this node or be from this input selection"""
        id: IntegerSocket
        """Secondary structure component to select"""

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
        id: InputInteger = 1,
    ):
        super().__init__(**{"And": and_, "Or": or_, "id": id})

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
        id = tree.inputs.integer(
            "id", 1, description="Secondary structure component to select"
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
                id, g.NamedAttribute.integer("sec_struct").o.attribute
            ),
        )

        group >> selection
        group.o.inverted >> inverted
