# Node group "MN Units" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import CustomGeometryGroup, FloatSocket, SocketAccessor
from nodebpy.types import InputFloat
from .mn_world_scale import MN_world_scale


class MNUnits(CustomGeometryGroup):
    """
    MN Units

    Parameters
    ----------
    value : InputFloat
        A value which will be scaled appropriately for the world

    Inputs
    ------
    i.value : FloatSocket
        A value which will be scaled appropriately for the world

    Outputs
    -------
    o.angstrom : FloatSocket
        Angstrom
    o.nanometre : FloatSocket
        Nanometre
    """

    _name = "MN Units"
    _color_tag = "CONVERTER"
    _tree_properties = {"node_tool_idname": "geometry.mn_units"}

    class _Inputs(SocketAccessor):
        value: FloatSocket
        """A value which will be scaled appropriately for the world"""

    class _Outputs(SocketAccessor):
        angstrom: FloatSocket
        """Angstrom"""
        nanometre: FloatSocket
        """Nanometre"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        value: InputFloat = 3.0,
    ):
        super().__init__(**{"Value": value})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        value = tree.inputs.float(
            "Value",
            3.0,
            description="A value which will be scaled appropriately for the world",
            min_value=-10_000.0,
            max_value=10_000.0,
        )
        angstrom = tree.outputs.float("Angstrom")
        nanometre = tree.outputs.float("Nanometre")

        math_1 = value * MN_world_scale()
        math_1 * 10.0 >> nanometre

        math_1 >> angstrom
