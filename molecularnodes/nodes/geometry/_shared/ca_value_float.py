# Node group "CA Value Float" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import CustomGeometryGroup, FloatSocket, SocketAccessor
from nodebpy.types import InputFloat
from ..atom_name import AtomName


class CAValueFloat(CustomGeometryGroup):
    """
    CA Value Float

    Parameters
    ----------
    value : InputFloat
        Value

    Inputs
    ------
    i.value : FloatSocket
        Value

    Outputs
    -------
    o.value : FloatSocket
        Value
    """

    _name = "CA Value Float"
    _color_tag = "CONVERTER"

    class _Inputs(SocketAccessor):
        value: FloatSocket
        """Value"""

    class _Outputs(SocketAccessor):
        value: FloatSocket
        """Value"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        value: InputFloat = 0.0,
    ):
        super().__init__(**{"Value": value})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        value = tree.inputs.float("Value", 0.0, hide_value=True)
        value_1 = tree.outputs.float("Value")

        g.IndexSwitch.float(AtomName(), (0.0, 0.0, value)) >> value_1
