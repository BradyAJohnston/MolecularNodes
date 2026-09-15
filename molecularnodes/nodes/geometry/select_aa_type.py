# Node-group asset "Select AA Type" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
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


class SelectAAType(AssetGeometryGroup):
    """
    Select AA Type

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
    o.is_polar : BooleanSocket
        The amino acid residue is polar
    o.is_apolar : BooleanSocket
        The amino acid residue is apolar
    o.is_acidic : BooleanSocket
        The amino acid residue is acidic
    o.is_basic : BooleanSocket
        The amino acid residue is basic
    o.is_aromatic : BooleanSocket
        The amino acid residue is aromatic
    """

    _name = "Select AA Type"
    _asset_name = "Select AA Type"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.select_aa_type"}

    class _Inputs(SocketAccessor):
        and_: BooleanSocket
        """The resulting selection must overlap with this input selection"""
        or_: BooleanSocket
        """The resulting selection can be calculated from this node or be from this input selection"""

    class _Outputs(SocketAccessor):
        is_polar: BooleanSocket
        """The amino acid residue is polar"""
        is_apolar: BooleanSocket
        """The amino acid residue is apolar"""
        is_acidic: BooleanSocket
        """The amino acid residue is acidic"""
        is_basic: BooleanSocket
        """The amino acid residue is basic"""
        is_aromatic: BooleanSocket
        """The amino acid residue is aromatic"""

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
        is_polar = tree.outputs.boolean(
            "is_polar", description="The amino acid residue is polar"
        )
        is_apolar = tree.outputs.boolean(
            "is_apolar", description="The amino acid residue is apolar"
        )
        is_acidic = tree.outputs.boolean(
            "is_acidic", description="The amino acid residue is acidic"
        )
        is_basic = tree.outputs.boolean(
            "is_basic", description="The amino acid residue is basic"
        )
        is_aromatic = tree.outputs.boolean(
            "is_aromatic", description="The amino acid residue is aromatic"
        )

        residue_name = ResidueName()
        with g.Frame("Residue Category"):
            index_switch = g.IndexSwitch.integer(
                residue_name,
                (
                    1,
                    2,
                    0,
                    3,
                    0,
                    3,
                    0,
                    1,
                    0,
                    1,
                    1,
                    2,
                    1,
                    4,
                    1,
                    0,
                    0,
                    4,
                    4,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                ),
            )
            _string = g.String(
                string="# ALA(0)\n# ARG(1)\n# ASN(2)\n# ASP(3)\n# CYS(4)\n# GLU(5)\n# GLN(6)\n# GLY(7)\n# HIS(8)\n# ILE(9)\n# LEU(10)\n# LYS(11)\n# MET(12)\n# PHE(13)\n# PRO(14)\n# SER(15)\n# THR(16)\n# TRP(17)\n# TYR(18)\n# VAL(19)"
            )
        switch = g.Compare.integer.equal(residue_name, -1).o.result.switch.integer(
            index_switch, -1
        )
        (
            BooleanAndOr(and_=and_, or_=or_, boolean=g.Compare.integer.equal(switch, 0))
            >> is_polar
        )
        (
            BooleanAndOr(and_=and_, or_=or_, boolean=g.Compare.integer.equal(switch, 1))
            >> is_apolar
        )
        (
            BooleanAndOr(and_=and_, or_=or_, boolean=g.Compare.integer.equal(switch, 3))
            >> is_acidic
        )
        (
            BooleanAndOr(and_=and_, or_=or_, boolean=g.Compare.integer.equal(switch, 2))
            >> is_basic
        )
        (
            BooleanAndOr(and_=and_, or_=or_, boolean=g.Compare.integer.equal(switch, 4))
            >> is_aromatic
        )


ASSET = SelectAAType

ASSET_METADATA = {
    "catalog_id": "bd1f205b-fea5-4700-b2c2-754f3321e969",
}
