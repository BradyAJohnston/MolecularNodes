# Node-group asset "Color pLDDT" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    AssetGeometryGroup,
    ColorSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputColor
from .b_factor import BFactor


class ColorPLDDT(AssetGeometryGroup):
    """
    Color pLDDT

    Parameters
    ----------
    _50 : InputColor
        Color for pLDTT < 50
    _70 : InputColor
        Color for 50 < pLDTT < 70
    socket_2 : InputColor
        Color for 70 < pLDTT < 90
    socket_3 : InputColor
        Color for 90 < pLDTT

    Inputs
    ------
    i._50 : ColorSocket
        Color for pLDTT < 50
    i._70 : ColorSocket
        Color for 50 < pLDTT < 70
    i.socket_2 : ColorSocket
        Color for 70 < pLDTT < 90
    i.socket_3 : ColorSocket
        Color for 90 < pLDTT

    Outputs
    -------
    o.color : ColorSocket
        Assigned color based on the pLDTT score
    """

    _name = "Color pLDDT"
    _asset_name = "Color pLDDT"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "COLOR"
    _tree_properties = {"node_tool_idname": "geometry.color_plddt"}

    class _Inputs(SocketAccessor):
        _50: ColorSocket
        """Color for pLDTT < 50"""
        _70: ColorSocket
        """Color for 50 < pLDTT < 70"""
        socket_2: ColorSocket
        """Color for 70 < pLDTT < 90"""
        socket_3: ColorSocket
        """Color for 90 < pLDTT"""

    class _Outputs(SocketAccessor):
        color: ColorSocket
        """Assigned color based on the pLDTT score"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        _50: InputColor = None,
        _70: InputColor = None,
        socket_2: InputColor = None,
        socket_3: InputColor = None,
    ):
        super().__init__(**{"<50": _50, "<70": _70, "<90": socket_2, ">90": socket_3})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        n_50 = tree.inputs.color(
            "<50",
            (1.000169, 0.20506974, 0.05950701, 1.0),
            description="Color for pLDTT < 50",
        )
        n_70 = tree.inputs.color(
            "<70",
            (1.0001687, 0.7083451, 0.006511816, 1.0),
            description="Color for 50 < pLDTT < 70",
        )
        n_90 = tree.inputs.color(
            "<90",
            (0.13015743, 0.5971759, 0.8962046, 1.0),
            description="Color for 70 < pLDTT < 90",
        )
        n_90_1 = tree.inputs.color(
            ">90", (0.0, 0.08649647, 0.6723945, 1.0), description="Color for 90 < pLDTT"
        )
        color = tree.outputs.color(
            "Color",
            (0.0, 0.0, 0.0, 0.0),
            description="Assigned color based on the pLDTT score",
        )

        group = BFactor()
        (
            (group > 90.0).switch.color(
                (group > 70.0).switch.color(
                    (group > 50.0).switch.color(n_50, n_70), n_90
                ),
                n_90_1,
            )
            >> color
        )


ASSET = ColorPLDDT

ASSET_METADATA = {
    "catalog_id": "d3f975df-8408-4972-a669-8187a57e01d0",
}
