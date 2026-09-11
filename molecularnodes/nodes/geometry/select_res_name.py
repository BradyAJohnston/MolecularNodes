# Node-group asset "Select Res Name" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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


class SelectResName(AssetGeometryGroup):
    """
    Select Res Name

    Parameters
    ----------
    and_ : InputBoolean
        The resulting selection must overlap with this input selection
    or_ : InputBoolean
        The resulting selection can be calculated from this node or be from this input selection
    ala : InputBoolean
        Select the residue ALA
    arg : InputBoolean
        Select the residue ARG
    asn : InputBoolean
        Select the residue ASN
    asp : InputBoolean
        Select the residue ASP
    cys : InputBoolean
        Select the residue CYS
    glu : InputBoolean
        Select the residue GLU
    gln : InputBoolean
        Select the residue GLN
    gly : InputBoolean
        Select the residue GLY
    his : InputBoolean
        Select the residue HIS
    ile : InputBoolean
        Select the residue ILE
    leu : InputBoolean
        Select the residue LEU
    lys : InputBoolean
        Select the residue LYS
    met : InputBoolean
        Select the residue MET
    phe : InputBoolean
        Select the residue PHE
    pro : InputBoolean
        Select the residue PRO
    ser : InputBoolean
        Select the residue SER
    thr : InputBoolean
        Select the residue THR
    trp : InputBoolean
        Select the residue TRP
    tyr : InputBoolean
        Select the residue TYR
    val : InputBoolean
        Select the residue VAL
    a : InputBoolean
        Select the residue A
    c : InputBoolean
        Select the residue C
    g : InputBoolean
        Select the residue G
    t : InputBoolean
        Select the residue T
    ra : InputBoolean
        Select the residue rA
    rc : InputBoolean
        Select the residue rC
    rg : InputBoolean
        Select the residue rG
    ru : InputBoolean
        Select the residue rU

    Inputs
    ------
    i.and_ : BooleanSocket
        The resulting selection must overlap with this input selection
    i.or_ : BooleanSocket
        The resulting selection can be calculated from this node or be from this input selection
    i.ala : BooleanSocket
        Select the residue ALA
    i.arg : BooleanSocket
        Select the residue ARG
    i.asn : BooleanSocket
        Select the residue ASN
    i.asp : BooleanSocket
        Select the residue ASP
    i.cys : BooleanSocket
        Select the residue CYS
    i.glu : BooleanSocket
        Select the residue GLU
    i.gln : BooleanSocket
        Select the residue GLN
    i.gly : BooleanSocket
        Select the residue GLY
    i.his : BooleanSocket
        Select the residue HIS
    i.ile : BooleanSocket
        Select the residue ILE
    i.leu : BooleanSocket
        Select the residue LEU
    i.lys : BooleanSocket
        Select the residue LYS
    i.met : BooleanSocket
        Select the residue MET
    i.phe : BooleanSocket
        Select the residue PHE
    i.pro : BooleanSocket
        Select the residue PRO
    i.ser : BooleanSocket
        Select the residue SER
    i.thr : BooleanSocket
        Select the residue THR
    i.trp : BooleanSocket
        Select the residue TRP
    i.tyr : BooleanSocket
        Select the residue TYR
    i.val : BooleanSocket
        Select the residue VAL
    i.a : BooleanSocket
        Select the residue A
    i.c : BooleanSocket
        Select the residue C
    i.g : BooleanSocket
        Select the residue G
    i.t : BooleanSocket
        Select the residue T
    i.ra : BooleanSocket
        Select the residue rA
    i.rc : BooleanSocket
        Select the residue rC
    i.rg : BooleanSocket
        Select the residue rG
    i.ru : BooleanSocket
        Select the residue rU

    Outputs
    -------
    o.selection : BooleanSocket
        The calculated selection
    o.inverted : BooleanSocket
        The inverse of the calculated selection
    """

    _name = "Select Res Name"
    _asset_name = "Select Res Name"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.select_res_name"}

    class _Inputs(SocketAccessor):
        and_: BooleanSocket
        """The resulting selection must overlap with this input selection"""
        or_: BooleanSocket
        """The resulting selection can be calculated from this node or be from this input selection"""
        ala: BooleanSocket
        """Select the residue ALA"""
        arg: BooleanSocket
        """Select the residue ARG"""
        asn: BooleanSocket
        """Select the residue ASN"""
        asp: BooleanSocket
        """Select the residue ASP"""
        cys: BooleanSocket
        """Select the residue CYS"""
        glu: BooleanSocket
        """Select the residue GLU"""
        gln: BooleanSocket
        """Select the residue GLN"""
        gly: BooleanSocket
        """Select the residue GLY"""
        his: BooleanSocket
        """Select the residue HIS"""
        ile: BooleanSocket
        """Select the residue ILE"""
        leu: BooleanSocket
        """Select the residue LEU"""
        lys: BooleanSocket
        """Select the residue LYS"""
        met: BooleanSocket
        """Select the residue MET"""
        phe: BooleanSocket
        """Select the residue PHE"""
        pro: BooleanSocket
        """Select the residue PRO"""
        ser: BooleanSocket
        """Select the residue SER"""
        thr: BooleanSocket
        """Select the residue THR"""
        trp: BooleanSocket
        """Select the residue TRP"""
        tyr: BooleanSocket
        """Select the residue TYR"""
        val: BooleanSocket
        """Select the residue VAL"""
        a: BooleanSocket
        """Select the residue A"""
        c: BooleanSocket
        """Select the residue C"""
        g: BooleanSocket
        """Select the residue G"""
        t: BooleanSocket
        """Select the residue T"""
        ra: BooleanSocket
        """Select the residue rA"""
        rc: BooleanSocket
        """Select the residue rC"""
        rg: BooleanSocket
        """Select the residue rG"""
        ru: BooleanSocket
        """Select the residue rU"""

    class _Outputs(SocketAccessor):
        selection: BooleanSocket
        """The calculated selection"""
        inverted: BooleanSocket
        """The inverse of the calculated selection"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        and_: InputBoolean = True,
        or_: InputBoolean = False,
        ala: InputBoolean = False,
        arg: InputBoolean = False,
        asn: InputBoolean = False,
        asp: InputBoolean = False,
        cys: InputBoolean = False,
        glu: InputBoolean = False,
        gln: InputBoolean = False,
        gly: InputBoolean = False,
        his: InputBoolean = False,
        ile: InputBoolean = False,
        leu: InputBoolean = False,
        lys: InputBoolean = False,
        met: InputBoolean = False,
        phe: InputBoolean = False,
        pro: InputBoolean = False,
        ser: InputBoolean = False,
        thr: InputBoolean = False,
        trp: InputBoolean = False,
        tyr: InputBoolean = False,
        val: InputBoolean = False,
        a: InputBoolean = False,
        c: InputBoolean = False,
        g: InputBoolean = False,
        t: InputBoolean = False,
        ra: InputBoolean = False,
        rc: InputBoolean = False,
        rg: InputBoolean = False,
        ru: InputBoolean = False,
    ):
        super().__init__(
            **{
                "And": and_,
                "Or": or_,
                "ALA": ala,
                "ARG": arg,
                "ASN": asn,
                "ASP": asp,
                "CYS": cys,
                "GLU": glu,
                "GLN": gln,
                "GLY": gly,
                "HIS": his,
                "ILE": ile,
                "LEU": leu,
                "LYS": lys,
                "MET": met,
                "PHE": phe,
                "PRO": pro,
                "SER": ser,
                "THR": thr,
                "TRP": trp,
                "TYR": tyr,
                "VAL": val,
                "A": a,
                "C": c,
                "G": g,
                "T": t,
                "rA": ra,
                "rC": rc,
                "rG": rg,
                "rU": ru,
            }
        )

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
        with tree.inputs.panel("Protein"):
            ala = tree.inputs.boolean(
                "ALA", False, description="Select the residue ALA"
            )
            arg = tree.inputs.boolean(
                "ARG", False, description="Select the residue ARG"
            )
            asn = tree.inputs.boolean(
                "ASN", False, description="Select the residue ASN"
            )
            asp = tree.inputs.boolean(
                "ASP", False, description="Select the residue ASP"
            )
            cys = tree.inputs.boolean(
                "CYS", False, description="Select the residue CYS"
            )
            glu = tree.inputs.boolean(
                "GLU", False, description="Select the residue GLU"
            )
            gln = tree.inputs.boolean(
                "GLN", False, description="Select the residue GLN"
            )
            gly = tree.inputs.boolean(
                "GLY", False, description="Select the residue GLY"
            )
            his = tree.inputs.boolean(
                "HIS", False, description="Select the residue HIS"
            )
            ile = tree.inputs.boolean(
                "ILE", False, description="Select the residue ILE"
            )
            leu = tree.inputs.boolean(
                "LEU", False, description="Select the residue LEU"
            )
            lys = tree.inputs.boolean(
                "LYS", False, description="Select the residue LYS"
            )
            met = tree.inputs.boolean(
                "MET", False, description="Select the residue MET"
            )
            phe = tree.inputs.boolean(
                "PHE", False, description="Select the residue PHE"
            )
            pro = tree.inputs.boolean(
                "PRO", False, description="Select the residue PRO"
            )
            ser = tree.inputs.boolean(
                "SER", False, description="Select the residue SER"
            )
            thr = tree.inputs.boolean(
                "THR", False, description="Select the residue THR"
            )
            trp = tree.inputs.boolean(
                "TRP", False, description="Select the residue TRP"
            )
            tyr = tree.inputs.boolean(
                "TYR", False, description="Select the residue TYR"
            )
            val = tree.inputs.boolean(
                "VAL", False, description="Select the residue VAL"
            )
        with tree.inputs.panel("DNA"):
            a = tree.inputs.boolean("A", False, description="Select the residue A")
            c_ = tree.inputs.boolean("C", False, description="Select the residue C")
            g_ = tree.inputs.boolean("G", False, description="Select the residue G")
            t = tree.inputs.boolean("T", False, description="Select the residue T")
        with tree.inputs.panel("RNA"):
            ra = tree.inputs.boolean("rA", False, description="Select the residue rA")
            rc = tree.inputs.boolean("rC", False, description="Select the residue rC")
            rg = tree.inputs.boolean("rG", False, description="Select the residue rG")
            ru = tree.inputs.boolean("rU", False, description="Select the residue rU")
        selection = tree.outputs.boolean(
            "Selection", description="The calculated selection"
        )
        inverted = tree.outputs.boolean(
            "Inverted", description="The inverse of the calculated selection"
        )

        index_switch = g.IndexSwitch.boolean(
            ResidueName(),
            (
                ala,
                arg,
                asn,
                asp,
                cys,
                glu,
                gln,
                gly,
                his,
                ile,
                leu,
                lys,
                met,
                phe,
                pro,
                ser,
                thr,
                trp,
                tyr,
                val,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                a,
                c_,
                g_,
                t,
                False,
                False,
                False,
                False,
                False,
                False,
                ra,
                rc,
                rg,
                ru,
            ),
        )
        group = BooleanAndOr(and_=and_, or_=or_, boolean=index_switch)

        group >> selection
        group.o.inverted >> inverted


ASSET = SelectResName

ASSET_METADATA = {
    "catalog_id": "bd1f205b-fea5-4700-b2c2-754f3321e969",
}
