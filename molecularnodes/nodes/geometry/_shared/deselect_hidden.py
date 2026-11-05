# Node group ".Deselect Hidden" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import BooleanSocket, CustomGeometryGroup, SocketAccessor
from nodebpy.types import InputBoolean


class DeselectHidden(CustomGeometryGroup):
    """
    .Deselect Hidden

    Parameters
    ----------
    selection : InputBoolean
        Selection of atoms for a style
    include_hidden : InputBoolean
        Keep atoms with the hidden attribute in the selection

    Inputs
    ------
    i.selection : BooleanSocket
        Selection of atoms for a style
    i.include_hidden : BooleanSocket
        Keep atoms with the hidden attribute in the selection

    Outputs
    -------
    o.selection : BooleanSocket
        The selection without hidden atoms. Atoms are not hidden when the attribute is missing
    """

    _name = ".Deselect Hidden"

    class _Inputs(SocketAccessor):
        selection: BooleanSocket
        """Selection of atoms for a style"""
        include_hidden: BooleanSocket
        """Keep atoms with the hidden attribute in the selection"""

    class _Outputs(SocketAccessor):
        selection: BooleanSocket
        """The selection without hidden atoms. Atoms are not hidden when the attribute is missing"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        selection: InputBoolean = True,
        include_hidden: InputBoolean = False,
    ):
        super().__init__(**{"Selection": selection, "Include Hidden": include_hidden})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Selection of atoms for a style",
            hide_value=True,
        )
        include_hidden = tree.inputs.boolean(
            "Include Hidden",
            False,
            description="Keep atoms with the hidden attribute in the selection",
        )
        selection_1 = tree.outputs.boolean(
            "Selection",
            description="The selection without hidden atoms. Atoms are not hidden when the attribute is missing",
        )

        (
            (
                selection
                & (include_hidden | ~g.NamedAttribute.boolean("hidden").o.attribute)
            )
            >> selection_1
        )
