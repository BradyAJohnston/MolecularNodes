# Node-group asset "Select Residue Class" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean
from .boolean_andor import BooleanAndOr
from .residue_name import ResidueName


class SelectResidueClass(AssetGeometryGroup):
    """
    Select Residue Class

    Parameters
    ----------
    and_ : InputBoolean
        The resulting selection must overlap with this input selection
    or_ : InputBoolean
        The resulting selection can be calculated from this node or be from this input selection

    Inputs
    ------
    i.and_ : BooleanSocket
        The resulting selection must overlap with this input selection
    i.or_ : BooleanSocket
        The resulting selection can be calculated from this node or be from this input selection

    Outputs
    -------
    o.positively_charged : BooleanSocket
        ARG, HIS and LYS
    o.negatively_charged : BooleanSocket
        ASP and GLU
    o.polar_uncharged : BooleanSocket
        SER, THR, ASN and GLN
    o.special_case : BooleanSocket
        CYS, SEC, GLY and PRO
    o.hydrophobic : BooleanSocket
        ALA, VAL, ILE, LEU, MET, PHE, TYR and TRP
    """

    _name = "Select Residue Class"
    _asset_name = "Select Residue Class"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.select_residue_class"}

    class _Inputs(SocketAccessor):
        and_: BooleanSocket
        """The resulting selection must overlap with this input selection"""
        or_: BooleanSocket
        """The resulting selection can be calculated from this node or be from this input selection"""

    class _Outputs(SocketAccessor):
        positively_charged: BooleanSocket
        """ARG, HIS and LYS"""
        negatively_charged: BooleanSocket
        """ASP and GLU"""
        polar_uncharged: BooleanSocket
        """SER, THR, ASN and GLN"""
        special_case: BooleanSocket
        """CYS, SEC, GLY and PRO"""
        hydrophobic: BooleanSocket
        """ALA, VAL, ILE, LEU, MET, PHE, TYR and TRP"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        and_: InputBoolean = True,
        or_: InputBoolean = False,
    ):
        super().__init__(**{"And": and_, "Or": or_})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        and_ = tree.inputs.boolean(
            "And",
            True,
            description="The resulting selection must overlap with this input selection",
            hide_value=True,
        )
        or_ = tree.inputs.boolean(
            "Or",
            False,
            description="The resulting selection can be calculated from this node or be from this input selection",
            hide_value=True,
        )
        positively_charged = tree.outputs.boolean(
            "Positively Charged", description="ARG, HIS and LYS"
        )
        negatively_charged = tree.outputs.boolean(
            "Negatively Charged", description="ASP and GLU"
        )
        polar_uncharged = tree.outputs.boolean(
            "Polar Uncharged", description="SER, THR, ASN and GLN"
        )
        special_case = tree.outputs.boolean(
            "Special Case", description="CYS, SEC, GLY and PRO"
        )
        hydrophobic = tree.outputs.boolean(
            "Hydrophobic", description="ALA, VAL, ILE, LEU, MET, PHE, TYR and TRP"
        )

        with g.Frame("Residue classes"):
            index_switch = g.IndexSwitch.integer(
                ResidueName(),
                (5, 1, 3, 2, 4, 2, 3, 4, 1, 5, 5, 1, 5, 5, 4, 3, 3, 5, 5, 5, 4),
            )
            _string = g.String(
                "Classes from the side chain chart on Wikipedia's Amino acid page: 1 positively charged, 2 negatively charged, 3 polar uncharged, 4 special case, 5 hydrophobic. Indexed by res_name, ALA (0) to VAL (19) then SEC (20); anything else falls outside the switch and gets 0."
            )
        (
            BooleanAndOr(
                and_=and_, or_=or_, boolean=g.Compare.integer.equal(index_switch, 1)
            )
            >> positively_charged
        )
        (
            BooleanAndOr(
                and_=and_, or_=or_, boolean=g.Compare.integer.equal(index_switch, 2)
            )
            >> negatively_charged
        )
        (
            BooleanAndOr(
                and_=and_, or_=or_, boolean=g.Compare.integer.equal(index_switch, 3)
            )
            >> polar_uncharged
        )
        (
            BooleanAndOr(
                and_=and_, or_=or_, boolean=g.Compare.integer.equal(index_switch, 4)
            )
            >> special_case
        )
        (
            BooleanAndOr(
                and_=and_, or_=or_, boolean=g.Compare.integer.equal(index_switch, 5)
            )
            >> hydrophobic
        )


ASSET = SelectResidueClass

ASSET_METADATA = {
    "description": "Select amino acids by residue class: positively charged, negatively charged, polar uncharged, special case or hydrophobic",
    "catalog_id": "bd1f205b-fea5-4700-b2c2-754f3321e969",
}
