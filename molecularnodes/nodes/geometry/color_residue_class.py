# Node-group asset "Color Residue Class" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
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
from .color import Color
from .select_residue_class import SelectResidueClass


class ColorResidueClass(AssetGeometryGroup):
    """
    Color Residue Class

    Parameters
    ----------
    positively_charged : InputColor
        Color for ARG, HIS and LYS
    negatively_charged : InputColor
        Color for ASP and GLU
    polar_uncharged : InputColor
        Color for SER, THR, ASN and GLN
    special_case : InputColor
        Color for CYS, SEC, GLY and PRO
    hydrophobic : InputColor
        Color for ALA, VAL, ILE, LEU, MET, PHE, TYR and TRP

    Inputs
    ------
    i.positively_charged : ColorSocket
        Color for ARG, HIS and LYS
    i.negatively_charged : ColorSocket
        Color for ASP and GLU
    i.polar_uncharged : ColorSocket
        Color for SER, THR, ASN and GLN
    i.special_case : ColorSocket
        Color for CYS, SEC, GLY and PRO
    i.hydrophobic : ColorSocket
        Color for ALA, VAL, ILE, LEU, MET, PHE, TYR and TRP

    Outputs
    -------
    o.color : ColorSocket
        Residue class colors, and the current color for anything that is not an amino acid
    """

    _name = "Color Residue Class"
    _asset_name = "Color Residue Class"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "COLOR"
    _tree_properties = {"node_tool_idname": "geometry.color_residue_class"}

    class _Inputs(SocketAccessor):
        positively_charged: ColorSocket
        """Color for ARG, HIS and LYS"""
        negatively_charged: ColorSocket
        """Color for ASP and GLU"""
        polar_uncharged: ColorSocket
        """Color for SER, THR, ASN and GLN"""
        special_case: ColorSocket
        """Color for CYS, SEC, GLY and PRO"""
        hydrophobic: ColorSocket
        """Color for ALA, VAL, ILE, LEU, MET, PHE, TYR and TRP"""

    class _Outputs(SocketAccessor):
        color: ColorSocket
        """Residue class colors, and the current color for anything that is not an amino acid"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        positively_charged: InputColor = None,
        negatively_charged: InputColor = None,
        polar_uncharged: InputColor = None,
        special_case: InputColor = None,
        hydrophobic: InputColor = None,
    ):
        super().__init__(
            **{
                "Positively Charged": positively_charged,
                "Negatively Charged": negatively_charged,
                "Polar Uncharged": polar_uncharged,
                "Special Case": special_case,
                "Hydrophobic": hydrophobic,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        positively_charged = tree.inputs.color(
            "Positively Charged",
            (0.13, 0.25, 0.85, 1.0),
            description="Color for ARG, HIS and LYS",
        )
        negatively_charged = tree.inputs.color(
            "Negatively Charged",
            (0.85, 0.08, 0.08, 1.0),
            description="Color for ASP and GLU",
        )
        polar_uncharged = tree.inputs.color(
            "Polar Uncharged",
            (0.25, 0.65, 0.85, 1.0),
            description="Color for SER, THR, ASN and GLN",
        )
        special_case = tree.inputs.color(
            "Special Case",
            (0.85, 0.65, 0.1, 1.0),
            description="Color for CYS, SEC, GLY and PRO",
        )
        hydrophobic = tree.inputs.color(
            "Hydrophobic",
            (0.8, 0.4, 0.15, 1.0),
            description="Color for ALA, VAL, ILE, LEU, MET, PHE, TYR and TRP",
        )
        color = tree.outputs.color(
            "Color",
            (1.0, 1.0, 1.0, 1.0),
            description="Residue class colors, and the current color for anything that is not an amino acid",
        )

        select_residue_class = SelectResidueClass()
        switch = select_residue_class.o.negatively_charged.switch.color(
            select_residue_class.o.positively_charged.switch.color(
                Color(), positively_charged
            ),
            negatively_charged,
        )
        switch_1 = select_residue_class.o.special_case.switch.color(
            select_residue_class.o.polar_uncharged.switch.color(
                switch, polar_uncharged
            ),
            special_case,
        )
        select_residue_class.o.hydrophobic.switch.color(switch_1, hydrophobic) >> color


ASSET = ColorResidueClass

ASSET_METADATA = {
    "description": "Color amino acids by residue class: positively charged, negatively charged, polar uncharged, special case or hydrophobic",
    "catalog_id": "d3f975df-8408-4972-a669-8187a57e01d0",
}
