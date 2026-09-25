# Node group "AngstromToWorld" (CompositorNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import CompositorNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import CustomCompositorGroup, FloatSocket, SocketAccessor
from nodebpy.types import InputFloat


class AngstromToWorld2(CustomCompositorGroup):
    """
    AngstromToWorld

    Parameters
    ----------
    angstrom : InputFloat
        Angstrom

    Inputs
    ------
    i.angstrom : FloatSocket
        Angstrom

    Outputs
    -------
    o.world : FloatSocket
        World
    """

    _name = "AngstromToWorld"

    class _Inputs(SocketAccessor):
        angstrom: FloatSocket
        """Angstrom"""

    class _Outputs(SocketAccessor):
        world: FloatSocket
        """World"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        angstrom: InputFloat = 0.5,
    ):
        super().__init__(**{"Angstrom": angstrom})

    def _build_group(self, tree: TreeBuilder[CompositorNodeTree]) -> None:
        angstrom = tree.inputs.float(
            "Angstrom", 0.5, min_value=-10_000.0, max_value=10_000.0
        )
        world = tree.outputs.float("World")

        angstrom * g.Value(0.1) >> world
