# Node group "Override Index" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    BooleanSocket,
    CustomGeometryGroup,
    IntegerSocket,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputInteger


class OverrideIndex(CustomGeometryGroup):
    """
    Override Index

    Parameters
    ----------
    selection : InputBoolean
        Selection
    index : InputInteger
        Index
    override : InputInteger
        Override

    Inputs
    ------
    i.selection : BooleanSocket
        Selection
    i.index : IntegerSocket
        Index
    i.override : IntegerSocket
        Override

    Outputs
    -------
    o.output : IntegerSocket
        Output
    """

    _name = "Override Index"
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        selection: BooleanSocket
        """Selection"""
        index: IntegerSocket
        """Index"""
        override: IntegerSocket
        """Override"""

    class _Outputs(SocketAccessor):
        output: IntegerSocket
        """Output"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        selection: InputBoolean = True,
        index: InputInteger = 0,
        override: InputInteger = 0,
    ):
        super().__init__(
            **{"Selection": selection, "Index": index, "Override": override}
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        selection = tree.inputs.boolean("Selection", True, hide_value=True)
        index = tree.inputs.integer("Index", 0, hide_value=True, default_input="INDEX")
        override = tree.inputs.integer("Override", 0)
        output = tree.outputs.integer("Output")

        selection.switch.integer(index, override) >> output
