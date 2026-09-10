# Node group '.MN_select_res_name_peptide' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import BooleanSocket, CustomGeometryGroup, SocketAccessor
from nodebpy.types import InputBoolean
from ..residue_name import ResidueName


class MN_select_res_name_peptide(CustomGeometryGroup):
    """
    .MN_select_res_name_peptide

    Parameters
    ----------
    ala : InputBoolean
        Select the AA residue ALA
    arg : InputBoolean
        Select the AA residue ARG
    asn : InputBoolean
        Select the AA residue ASN
    asp : InputBoolean
        Select the AA residue ASP
    cys : InputBoolean
        Select the AA residue CYS
    glu : InputBoolean
        Select the AA residue GLU
    gln : InputBoolean
        Select the AA residue GLN
    gly : InputBoolean
        Select the AA residue GLY
    his : InputBoolean
        Select the AA residue HIS
    ile : InputBoolean
        Select the AA residue ILE
    leu : InputBoolean
        Select the AA residue LEU
    lys : InputBoolean
        Select the AA residue LYS
    met : InputBoolean
        Select the AA residue MET
    phe : InputBoolean
        Select the AA residue PHE
    pro : InputBoolean
        Select the AA residue PRO
    ser : InputBoolean
        Select the AA residue SER
    thr : InputBoolean
        Select the AA residue THR
    trp : InputBoolean
        Select the AA residue TRP
    tyr : InputBoolean
        Select the AA residue TYR
    val : InputBoolean
        Select the AA residue VAL

    Inputs
    ------
    i.ala : BooleanSocket
        Select the AA residue ALA
    i.arg : BooleanSocket
        Select the AA residue ARG
    i.asn : BooleanSocket
        Select the AA residue ASN
    i.asp : BooleanSocket
        Select the AA residue ASP
    i.cys : BooleanSocket
        Select the AA residue CYS
    i.glu : BooleanSocket
        Select the AA residue GLU
    i.gln : BooleanSocket
        Select the AA residue GLN
    i.gly : BooleanSocket
        Select the AA residue GLY
    i.his : BooleanSocket
        Select the AA residue HIS
    i.ile : BooleanSocket
        Select the AA residue ILE
    i.leu : BooleanSocket
        Select the AA residue LEU
    i.lys : BooleanSocket
        Select the AA residue LYS
    i.met : BooleanSocket
        Select the AA residue MET
    i.phe : BooleanSocket
        Select the AA residue PHE
    i.pro : BooleanSocket
        Select the AA residue PRO
    i.ser : BooleanSocket
        Select the AA residue SER
    i.thr : BooleanSocket
        Select the AA residue THR
    i.trp : BooleanSocket
        Select the AA residue TRP
    i.tyr : BooleanSocket
        Select the AA residue TYR
    i.val : BooleanSocket
        Select the AA residue VAL

    Outputs
    -------
    o.selection : BooleanSocket
        The calculated selection
    o.inverted : BooleanSocket
        The inverse of the calculated selection
    """

    _name = ".MN_select_res_name_peptide"
    _tree_properties = {"node_tool_idname": "geometry._mn_select_res_name_peptide"}

    class _Inputs(SocketAccessor):
        ala: BooleanSocket
        """Select the AA residue ALA"""
        arg: BooleanSocket
        """Select the AA residue ARG"""
        asn: BooleanSocket
        """Select the AA residue ASN"""
        asp: BooleanSocket
        """Select the AA residue ASP"""
        cys: BooleanSocket
        """Select the AA residue CYS"""
        glu: BooleanSocket
        """Select the AA residue GLU"""
        gln: BooleanSocket
        """Select the AA residue GLN"""
        gly: BooleanSocket
        """Select the AA residue GLY"""
        his: BooleanSocket
        """Select the AA residue HIS"""
        ile: BooleanSocket
        """Select the AA residue ILE"""
        leu: BooleanSocket
        """Select the AA residue LEU"""
        lys: BooleanSocket
        """Select the AA residue LYS"""
        met: BooleanSocket
        """Select the AA residue MET"""
        phe: BooleanSocket
        """Select the AA residue PHE"""
        pro: BooleanSocket
        """Select the AA residue PRO"""
        ser: BooleanSocket
        """Select the AA residue SER"""
        thr: BooleanSocket
        """Select the AA residue THR"""
        trp: BooleanSocket
        """Select the AA residue TRP"""
        tyr: BooleanSocket
        """Select the AA residue TYR"""
        val: BooleanSocket
        """Select the AA residue VAL"""

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
    ):
        super().__init__(
            **{
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
            }
        )

    def _build_group(self, tree):
        ala = tree.inputs.boolean("ALA", False, description="Select the AA residue ALA")
        arg = tree.inputs.boolean("ARG", False, description="Select the AA residue ARG")
        asn = tree.inputs.boolean("ASN", False, description="Select the AA residue ASN")
        asp = tree.inputs.boolean("ASP", False, description="Select the AA residue ASP")
        cys = tree.inputs.boolean("CYS", False, description="Select the AA residue CYS")
        glu = tree.inputs.boolean("GLU", False, description="Select the AA residue GLU")
        gln = tree.inputs.boolean("GLN", False, description="Select the AA residue GLN")
        gly = tree.inputs.boolean("GLY", False, description="Select the AA residue GLY")
        his = tree.inputs.boolean("HIS", False, description="Select the AA residue HIS")
        ile = tree.inputs.boolean("ILE", False, description="Select the AA residue ILE")
        leu = tree.inputs.boolean("LEU", False, description="Select the AA residue LEU")
        lys = tree.inputs.boolean("LYS", False, description="Select the AA residue LYS")
        met = tree.inputs.boolean("MET", False, description="Select the AA residue MET")
        phe = tree.inputs.boolean("PHE", False, description="Select the AA residue PHE")
        pro = tree.inputs.boolean("PRO", False, description="Select the AA residue PRO")
        ser = tree.inputs.boolean("SER", False, description="Select the AA residue SER")
        thr = tree.inputs.boolean("THR", False, description="Select the AA residue THR")
        trp = tree.inputs.boolean("TRP", False, description="Select the AA residue TRP")
        tyr = tree.inputs.boolean("TYR", False, description="Select the AA residue TYR")
        val = tree.inputs.boolean("VAL", False, description="Select the AA residue VAL")
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
            ),
        )
        ~index_switch.o.output >> inverted

        index_switch >> selection
