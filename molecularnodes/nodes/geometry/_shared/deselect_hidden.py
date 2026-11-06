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

    Inputs
    ------
    i.selection : BooleanSocket
        Selection of atoms for a style

    Outputs
    -------
    o.selection : BooleanSocket
        The selection without hidden atoms. Atoms are not hidden when the attribute is missing
    """

    _name = ".Deselect Hidden"

    class _Inputs(SocketAccessor):
        selection: BooleanSocket
        """Selection of atoms for a style"""

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
    ):
        super().__init__(**{"Selection": selection})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Selection of atoms for a style",
            hide_value=True,
        )
        selection_1 = tree.outputs.boolean(
            "Selection",
            description="The selection without hidden atoms. Atoms are not hidden when the attribute is missing",
        )

        (selection & ~g.NamedAttribute.boolean("hidden").o.attribute) >> selection_1
