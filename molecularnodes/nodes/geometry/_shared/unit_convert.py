# Node group 'Unit Convert' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from nodebpy import geometry as g
from nodebpy.builder import CustomGeometryGroup, FloatSocket, MenuSocket, SocketAccessor
from nodebpy.types import InputFloat, InputMenu
from .mn_world_scale import MN_world_scale


class UnitConvert(CustomGeometryGroup):
    """
    Unit Convert

    Parameters
    ----------
    distance_type : InputMenu | Literal["Angstrom", "Nanometre", "Micrometre"]
        What unit to scale the value to
    from_ : InputFloat
        A value which will be scaled appropriately for the world

    Inputs
    ------
    i.distance_type : MenuSocket
        What unit to scale the value to
    i.from_ : FloatSocket
        A value which will be scaled appropriately for the world

    Outputs
    -------
    o.world : FloatSocket
        The value that has been scaled appropriately for the world space
    """

    _name = "Unit Convert"
    _color_tag = "CONVERTER"
    _tree_properties = {"node_tool_idname": "geometry.unit_convert"}

    class _Inputs(SocketAccessor):
        distance_type: MenuSocket
        """What unit to scale the value to"""
        from_: FloatSocket
        """A value which will be scaled appropriately for the world"""

    class _Outputs(SocketAccessor):
        world: FloatSocket
        """The value that has been scaled appropriately for the world space"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        distance_type: InputMenu
        | Literal["Angstrom", "Nanometre", "Micrometre"] = "Angstrom",
        from_: InputFloat = 3.0,
    ):
        super().__init__(**{"Distance Type": distance_type, "From": from_})

    def _build_group(self, tree):
        distance_type = tree.inputs.menu(
            "Distance Type",
            description="What unit to scale the value to",
            optional_label=True,
        )
        from_ = tree.inputs.float(
            "From",
            3.0,
            description="A value which will be scaled appropriately for the world",
            min_value=-10_000.0,
            max_value=10_000.0,
        )
        world = tree.outputs.float(
            "World",
            description="The value that has been scaled appropriately for the world space",
        )

        math_1 = from_ * MN_world_scale()
        math_2 = math_1 * 10.0
        (
            g.MenuSwitch.float(
                distance_type,
                {
                    "Angstrom": math_1,
                    "Nanometre": math_2,
                    "Micrometre": math_2 * 1000.0,
                },
            )
            >> world
        )

        distance_type.default_value = "Angstrom"
