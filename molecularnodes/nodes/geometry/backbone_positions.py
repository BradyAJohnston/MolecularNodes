# Node-group asset "Backbone Positions" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputMenu
from .backbone_c import BackboneC
from .backbone_ca import BackboneCA
from .backbone_n import BackboneN
from .backbone_nh import BackboneNH
from .backbone_o import BackboneO


class BackbonePositions(AssetGeometryGroup):
    """
    Backbone Positions

    Parameters
    ----------
    method : InputMenu | Literal["Read", "Compute"]
        Method

    Inputs
    ------
    i.method : MenuSocket
        Method

    Outputs
    -------
    o.o : VectorSocket
        The position of the backbone _O_ atom for the residue
    o.c : VectorSocket
        The position of the backbone _C_ atom for the residue
    o.ca : VectorSocket
        The position of the backbone _CA_ atom for the residue
    o.n : VectorSocket
        The position of the backbone _N_ atom for the residue
    o.nh : VectorSocket
        NH
    """

    _name = "Backbone Positions"
    _asset_name = "Backbone Positions"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.backbone_positions"}

    class _Inputs(SocketAccessor):
        method: MenuSocket
        """Method"""

    class _Outputs(SocketAccessor):
        o: VectorSocket
        """The position of the backbone _O_ atom for the residue"""
        c: VectorSocket
        """The position of the backbone _C_ atom for the residue"""
        ca: VectorSocket
        """The position of the backbone _CA_ atom for the residue"""
        n: VectorSocket
        """The position of the backbone _N_ atom for the residue"""
        nh: VectorSocket
        """NH"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        method: InputMenu | Literal["Read", "Compute"] = "Compute",
    ):
        super().__init__(**{"Method": method})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        method = tree.inputs.menu("Method", expanded=True, optional_label=True)
        o = tree.outputs.vector(
            "O", description="The position of the backbone _O_ atom for the residue"
        )
        c_ = tree.outputs.vector(
            "C", description="The position of the backbone _C_ atom for the residue"
        )
        ca = tree.outputs.vector(
            "CA", description="The position of the backbone _CA_ atom for the residue"
        )
        n = tree.outputs.vector(
            "N", description="The position of the backbone _N_ atom for the residue"
        )
        nh = tree.outputs.vector("NH")

        menu_switch = g.MenuSwitch.integer(method, {"Read": 0, "Compute": 1})
        (
            BackboneO(
                method=g.IndexSwitch.menu(menu_switch.o.output, ("Read", "Compute"))
            )
            >> o
        )
        (
            BackboneC(
                method=g.IndexSwitch.menu(menu_switch.o.output, ("Read", "Compute"))
            )
            >> c_
        )
        (
            BackboneCA(
                method=g.IndexSwitch.menu(menu_switch.o.output, ("Read", "Compute"))
            )
            >> ca
        )
        (
            BackboneN(
                method=g.IndexSwitch.menu(menu_switch.o.output, ("Read", "Compute"))
            )
            >> n
        )
        (
            BackboneNH(
                menu=g.IndexSwitch.menu(menu_switch.o.output, ("Read", "Compute"))
            )
            >> nh
        )

        method.default_value = "Compute"


ASSET = BackbonePositions

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
