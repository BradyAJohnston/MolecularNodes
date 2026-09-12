# Node group ".MN_select_sec_struct" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import BooleanSocket, CustomGeometryGroup, SocketAccessor
from nodebpy.types import InputBoolean
from ..is_helix import IsHelix
from ..is_loop import IsLoop
from ..is_sheet import IsSheet


class MN_select_sec_struct(CustomGeometryGroup):
    """
    .MN_select_sec_struct

    Parameters
    ----------
    and_ : InputBoolean
        The resulting selection must overlap with this input selection

    Inputs
    ------
    i.and_ : BooleanSocket
        The resulting selection must overlap with this input selection

    Outputs
    -------
    o.is_helix : BooleanSocket
        Is Helix
    o.is_sheet : BooleanSocket
        Is Sheet
    o.is_structured : BooleanSocket
        Is Structured
    o.is_loop : BooleanSocket
        Is Loop
    """

    _name = ".MN_select_sec_struct"
    _tree_properties = {"node_tool_idname": "geometry._mn_select_sec_struct"}

    class _Inputs(SocketAccessor):
        and_: BooleanSocket
        """The resulting selection must overlap with this input selection"""

    class _Outputs(SocketAccessor):
        is_helix: BooleanSocket
        """Is Helix"""
        is_sheet: BooleanSocket
        """Is Sheet"""
        is_structured: BooleanSocket
        """Is Structured"""
        is_loop: BooleanSocket
        """Is Loop"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        and_: InputBoolean = True,
    ):
        super().__init__(**{"And": and_})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        and_ = tree.inputs.boolean(
            "And",
            True,
            description="The resulting selection must overlap with this input selection",
            hide_value=True,
        )
        is_helix = tree.outputs.boolean("Is Helix")
        is_sheet = tree.outputs.boolean("Is Sheet")
        is_structured = tree.outputs.boolean("Is Structured")
        is_loop = tree.outputs.boolean("Is Loop")

        IsHelix(and_=and_) >> is_helix
        group = IsLoop(and_=and_)
        ~group.o.selection >> is_structured
        IsSheet(and_=and_) >> is_sheet

        group >> is_loop
