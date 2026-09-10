# Node-group asset 'Fractionate Float' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    IntegerSocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputFloat, InputMenu


class FractionateFloat(AssetGeometryGroup):
    """
    Fractionate Float

    Parameters
    ----------
    menu : InputMenu | Literal["Linear", "Smoother"]
        Menu
    value : InputFloat
        The value to fractionate

    Inputs
    ------
    i.menu : MenuSocket
        Menu
    i.value : FloatSocket
        The value to fractionate

    Outputs
    -------
    o.fraction : FloatSocket
        Fractional component of the value, between 0 and 1
    o.floor : IntegerSocket
        The floor of the value; the integer rounded down
    o.ceiling : IntegerSocket
        The ceiling of the value, the integer rounded up
    """

    _name = "Fractionate Float"
    _asset_name = "Fractionate Float"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {"node_tool_idname": "geometry.fractionate_float"}

    class _Inputs(SocketAccessor):
        menu: MenuSocket
        """Menu"""
        value: FloatSocket
        """The value to fractionate"""

    class _Outputs(SocketAccessor):
        fraction: FloatSocket
        """Fractional component of the value, between 0 and 1"""
        floor: IntegerSocket
        """The floor of the value; the integer rounded down"""
        ceiling: IntegerSocket
        """The ceiling of the value, the integer rounded up"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        menu: InputMenu | Literal["Linear", "Smoother"] = "Linear",
        value: InputFloat = 0.0,
    ):
        super().__init__(**{"Menu": menu, "Value": value})

    def _build_group(self, tree):
        menu = tree.inputs.menu("Menu", expanded=True, optional_label=True)
        value = tree.inputs.float("Value", 0.0, description="The value to fractionate")
        fraction = tree.outputs.float(
            "Fraction", description="Fractional component of the value, between 0 and 1"
        )
        floor = tree.outputs.integer(
            "Floor", description="The floor of the value; the integer rounded down"
        )
        ceiling = tree.outputs.integer(
            "Ceiling", description="The ceiling of the value, the integer rounded up"
        )

        float_to_integer = g.FloatToInteger(float=value, rounding_mode="FLOOR")
        float_to_integer_1 = g.FloatToInteger(float=value, rounding_mode="CEILING")
        math_1 = value.fraction()
        (
            g.MenuSwitch.float(
                menu,
                {
                    "Linear": math_1,
                    "Smoother": math_1.map_range(interpolation_type="SMOOTHERSTEP"),
                },
            )
            >> fraction
        )

        float_to_integer >> floor
        float_to_integer_1 >> ceiling

        menu.default_value = "Linear"


ASSET = FractionateFloat

ASSET_METADATA = {
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
