# Node-group asset 'Color Common' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    ColorSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputColor
from .atomic_number import AtomicNumber
from .between_integer import BetweenInteger
from .color import Color


class ColorCommon(AssetGeometryGroup):
    """
    Color Common

    Parameters
    ----------
    hydrogen : InputColor
        Color to set for the element Hydrogen
    carbon : InputColor
        Color to set for the element Carbon
    nitrogen : InputColor
        Color to set for the element Nitrogen
    oxygen : InputColor
        Color to set for the element Oxygen
    phosphorous : InputColor
        Color to set for the element Phosphorous
    sulfur : InputColor
        Color to set for the element Sulfur

    Inputs
    ------
    i.hydrogen : ColorSocket
        Color to set for the element Hydrogen
    i.carbon : ColorSocket
        Color to set for the element Carbon
    i.nitrogen : ColorSocket
        Color to set for the element Nitrogen
    i.oxygen : ColorSocket
        Color to set for the element Oxygen
    i.phosphorous : ColorSocket
        Color to set for the element Phosphorous
    i.sulfur : ColorSocket
        Color to set for the element Sulfur

    Outputs
    -------
    o.color : ColorSocket
        The output colors for the common elements
    """

    _name = "Color Common"
    _asset_name = "Color Common"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "COLOR"
    _tree_properties = {"node_tool_idname": "geometry.color_common"}

    class _Inputs(SocketAccessor):
        hydrogen: ColorSocket
        """Color to set for the element Hydrogen"""
        carbon: ColorSocket
        """Color to set for the element Carbon"""
        nitrogen: ColorSocket
        """Color to set for the element Nitrogen"""
        oxygen: ColorSocket
        """Color to set for the element Oxygen"""
        phosphorous: ColorSocket
        """Color to set for the element Phosphorous"""
        sulfur: ColorSocket
        """Color to set for the element Sulfur"""

    class _Outputs(SocketAccessor):
        color: ColorSocket
        """The output colors for the common elements"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        hydrogen: InputColor = None,
        carbon: InputColor = None,
        nitrogen: InputColor = None,
        oxygen: InputColor = None,
        phosphorous: InputColor = None,
        sulfur: InputColor = None,
    ):
        super().__init__(
            **{
                "Hydrogen": hydrogen,
                "Carbon": carbon,
                "Nitrogen": nitrogen,
                "Oxygen": oxygen,
                "Phosphorous": phosphorous,
                "Sulfur": sulfur,
            }
        )

    def _build_group(self, tree):
        hydrogen = tree.inputs.color(
            "Hydrogen",
            (1.0, 1.0, 1.0, 1.0),
            description="Color to set for the element Hydrogen",
        )
        carbon = tree.inputs.color(
            "Carbon",
            (0.20190106, 0.20190106, 0.20190106, 1.0),
            description="Color to set for the element Carbon",
        )
        nitrogen = tree.inputs.color(
            "Nitrogen",
            (0.16, 0.23333497, 0.8, 1.0),
            description="Color to set for the element Nitrogen",
        )
        oxygen = tree.inputs.color(
            "Oxygen",
            (0.8, 0.1610207, 0.16, 1.0),
            description="Color to set for the element Oxygen",
        )
        phosphorous = tree.inputs.color(
            "Phosphorous",
            (0.8, 0.17181273, 0.5252497, 1.0),
            description="Color to set for the element Phosphorous",
        )
        sulfur = tree.inputs.color(
            "Sulfur",
            (0.8, 0.722058, 0.05199071, 1.0),
            description="Color to set for the element Sulfur",
        )
        color = tree.outputs.color(
            "Color",
            (0.0, 0.0, 0.0, 0.0),
            description="The output colors for the common elements",
        )

        group = Color()
        group_1 = AtomicNumber()
        index_switch = g.IndexSwitch.color(
            group_1,
            (
                group,
                hydrogen,
                group,
                group,
                group,
                group,
                carbon,
                nitrogen,
                oxygen,
                group,
                group,
                group,
                group,
                group,
                group,
                phosphorous,
                sulfur,
            ),
        )
        (
            BetweenInteger(value=group_1, upper=16).o.boolean.switch.color(
                group, index_switch
            )
            >> color
        )


ASSET = ColorCommon

ASSET_METADATA = {
    "catalog_id": "d3f975df-8408-4972-a669-8187a57e01d0",
}
