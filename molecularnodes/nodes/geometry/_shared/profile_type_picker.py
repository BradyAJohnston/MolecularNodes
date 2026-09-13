# Node group ".Profile Type Picker" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    CustomGeometryGroup,
    IntegerSocket,
    MenuSocket,
    SocketAccessor,
)
from nodebpy.types import InputMenu


class ProfileTypePicker(CustomGeometryGroup):
    """
    .Profile Type Picker

    Parameters
    ----------
    menu : InputMenu | Literal["Default Profile", "Custom Profile"]
        Menu

    Inputs
    ------
    i.menu : MenuSocket
        Menu

    Outputs
    -------
    o.output : IntegerSocket
        Output
    """

    _name = ".Profile Type Picker"
    _tree_properties = {"node_tool_idname": "geometry._profile_type_picker"}

    class _Inputs(SocketAccessor):
        menu: MenuSocket
        """Menu"""

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
        menu: InputMenu
        | Literal["Default Profile", "Custom Profile"] = "Custom Profile",
    ):
        super().__init__(**{"Menu": menu})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        menu = tree.inputs.menu("Menu", optional_label=True)
        output = tree.outputs.integer("Output")

        (
            g.MenuSwitch.integer(menu, {"Default Profile": 0, "Custom Profile": 1})
            >> output
        )

        menu.default_value = "Custom Profile"
