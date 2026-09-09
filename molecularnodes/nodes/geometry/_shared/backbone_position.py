# Node group '.Backbone Position' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from nodebpy import geometry as g
from nodebpy.builder import (
    CustomGeometryGroup,
    MenuSocket,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputMenu
from ..fallback_vector import FallbackVector
from ..residue_mask import ResidueMask


class BackbonePosition(CustomGeometryGroup):
    """
    .Backbone Position

    Parameters
    ----------
    method : InputMenu | Literal["Read", "Compute"]
        Method
    menu : InputMenu | Literal["backbone_N", "backbone_CA", "backbone_C", "backbone_O"]
        the particular backbone residue to read the value from

    Inputs
    ------
    i.method : MenuSocket
        Method
    i.menu : MenuSocket
        the particular backbone residue to read the value from

    Outputs
    -------
    o.position : VectorSocket
        Position
    """

    _name = ".Backbone Position"
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry._backbone_position"}

    class _Inputs(SocketAccessor):
        method: MenuSocket
        """Method"""
        menu: MenuSocket
        """the particular backbone residue to read the value from"""

    class _Outputs(SocketAccessor):
        position: VectorSocket
        """Position"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        method: InputMenu | Literal["Read", "Compute"] = "Read",
        menu: InputMenu
        | Literal[
            "backbone_N", "backbone_CA", "backbone_C", "backbone_O"
        ] = "backbone_N",
    ):
        super().__init__(**{"Method": method, "Menu": menu})

    def _build_group(self, tree):
        method = tree.inputs.menu("Method", expanded=True, optional_label=True)
        menu = tree.inputs.menu(
            "Menu",
            description="the particular backbone residue to read the value from",
            optional_label=True,
        )
        position = tree.outputs.vector("Position")

        menu_switch = g.MenuSwitch.integer(
            menu, {"backbone_N": 1, "backbone_CA": 2, "backbone_C": 3, "backbone_O": 4}
        )
        index_switch = g.IndexSwitch.string(
            menu_switch.o.output,
            ("", "backbone_N", "backbone_CA", "backbone_C", "backbone_O"),
        )
        group = FallbackVector(
            name=index_switch,
            fallback=ResidueMask(atom_name=menu_switch.o.output).o.position,
        )
        (
            g.MenuSwitch.vector(
                method,
                {
                    "Read": g.NamedAttribute.vector(index_switch).o.attribute,
                    "Compute": group.o.output.point.at(
                        ResidueMask(atom_name=2).o.index
                    ),
                },
            )
            >> position
        )

        method.default_value = "Read"
        menu.default_value = "backbone_N"
