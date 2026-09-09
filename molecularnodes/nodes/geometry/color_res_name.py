# Node-group asset 'Color Res Name' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from .between_integer import BetweenInteger
from .color import Color
from .residue_name import ResidueName


class ColorResName(AssetGeometryGroup):
    """
    Color Res Name

    Parameters
    ----------
    ala : InputColor
        Set the color for the residue ALA
    arg : InputColor
        Set the color for the residue ARG
    asn : InputColor
        Set the color for the residue ASN
    asp : InputColor
        Set the color for the residue ASP
    cys : InputColor
        Set the color for the residue CYS
    glu : InputColor
        Set the color for the residue GLU
    gln : InputColor
        Set the color for the residue GLN
    gly : InputColor
        Set the color for the residue GLY
    his : InputColor
        Set the color for the residue HIS
    ile : InputColor
        Set the color for the residue ILE
    leu : InputColor
        Set the color for the residue LEU
    lys : InputColor
        Set the color for the residue LYS
    met : InputColor
        Set the color for the residue MET
    phe : InputColor
        Set the color for the residue PHE
    pro : InputColor
        Set the color for the residue PRO
    ser : InputColor
        Set the color for the residue SER
    thr : InputColor
        Set the color for the residue THR
    trp : InputColor
        Set the color for the residue TRP
    tyr : InputColor
        Set the color for the residue TYR
    val : InputColor
        Set the color for the residue VAL
    a : InputColor
        Set the color for the residue A
    c : InputColor
        Set the color for the residue C
    g : InputColor
        Set the color for the residue G
    t : InputColor
        Set the color for the residue T
    ra : InputColor
        Set the color for the residue rA
    rc : InputColor
        Set the color for the residue rC
    rg : InputColor
        Set the color for the residue rG
    ru : InputColor
        Set the color for the residue rU

    Inputs
    ------
    i.ala : ColorSocket
        Set the color for the residue ALA
    i.arg : ColorSocket
        Set the color for the residue ARG
    i.asn : ColorSocket
        Set the color for the residue ASN
    i.asp : ColorSocket
        Set the color for the residue ASP
    i.cys : ColorSocket
        Set the color for the residue CYS
    i.glu : ColorSocket
        Set the color for the residue GLU
    i.gln : ColorSocket
        Set the color for the residue GLN
    i.gly : ColorSocket
        Set the color for the residue GLY
    i.his : ColorSocket
        Set the color for the residue HIS
    i.ile : ColorSocket
        Set the color for the residue ILE
    i.leu : ColorSocket
        Set the color for the residue LEU
    i.lys : ColorSocket
        Set the color for the residue LYS
    i.met : ColorSocket
        Set the color for the residue MET
    i.phe : ColorSocket
        Set the color for the residue PHE
    i.pro : ColorSocket
        Set the color for the residue PRO
    i.ser : ColorSocket
        Set the color for the residue SER
    i.thr : ColorSocket
        Set the color for the residue THR
    i.trp : ColorSocket
        Set the color for the residue TRP
    i.tyr : ColorSocket
        Set the color for the residue TYR
    i.val : ColorSocket
        Set the color for the residue VAL
    i.a : ColorSocket
        Set the color for the residue A
    i.c : ColorSocket
        Set the color for the residue C
    i.g : ColorSocket
        Set the color for the residue G
    i.t : ColorSocket
        Set the color for the residue T
    i.ra : ColorSocket
        Set the color for the residue rA
    i.rc : ColorSocket
        Set the color for the residue rC
    i.rg : ColorSocket
        Set the color for the residue rG
    i.ru : ColorSocket
        Set the color for the residue rU

    Outputs
    -------
    o.color : ColorSocket
        Color
    """

    _name = "Color Res Name"
    _asset_name = "Color Res Name"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "COLOR"
    _tree_properties = {"node_tool_idname": "geometry.color_res_name"}

    class _Inputs(SocketAccessor):
        ala: ColorSocket
        """Set the color for the residue ALA"""
        arg: ColorSocket
        """Set the color for the residue ARG"""
        asn: ColorSocket
        """Set the color for the residue ASN"""
        asp: ColorSocket
        """Set the color for the residue ASP"""
        cys: ColorSocket
        """Set the color for the residue CYS"""
        glu: ColorSocket
        """Set the color for the residue GLU"""
        gln: ColorSocket
        """Set the color for the residue GLN"""
        gly: ColorSocket
        """Set the color for the residue GLY"""
        his: ColorSocket
        """Set the color for the residue HIS"""
        ile: ColorSocket
        """Set the color for the residue ILE"""
        leu: ColorSocket
        """Set the color for the residue LEU"""
        lys: ColorSocket
        """Set the color for the residue LYS"""
        met: ColorSocket
        """Set the color for the residue MET"""
        phe: ColorSocket
        """Set the color for the residue PHE"""
        pro: ColorSocket
        """Set the color for the residue PRO"""
        ser: ColorSocket
        """Set the color for the residue SER"""
        thr: ColorSocket
        """Set the color for the residue THR"""
        trp: ColorSocket
        """Set the color for the residue TRP"""
        tyr: ColorSocket
        """Set the color for the residue TYR"""
        val: ColorSocket
        """Set the color for the residue VAL"""
        a: ColorSocket
        """Set the color for the residue A"""
        c: ColorSocket
        """Set the color for the residue C"""
        g: ColorSocket
        """Set the color for the residue G"""
        t: ColorSocket
        """Set the color for the residue T"""
        ra: ColorSocket
        """Set the color for the residue rA"""
        rc: ColorSocket
        """Set the color for the residue rC"""
        rg: ColorSocket
        """Set the color for the residue rG"""
        ru: ColorSocket
        """Set the color for the residue rU"""

    class _Outputs(SocketAccessor):
        color: ColorSocket
        """Color"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        ala: InputColor = None,
        arg: InputColor = None,
        asn: InputColor = None,
        asp: InputColor = None,
        cys: InputColor = None,
        glu: InputColor = None,
        gln: InputColor = None,
        gly: InputColor = None,
        his: InputColor = None,
        ile: InputColor = None,
        leu: InputColor = None,
        lys: InputColor = None,
        met: InputColor = None,
        phe: InputColor = None,
        pro: InputColor = None,
        ser: InputColor = None,
        thr: InputColor = None,
        trp: InputColor = None,
        tyr: InputColor = None,
        val: InputColor = None,
        a: InputColor = None,
        c: InputColor = None,
        g: InputColor = None,
        t: InputColor = None,
        ra: InputColor = None,
        rc: InputColor = None,
        rg: InputColor = None,
        ru: InputColor = None,
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

    def _build_group(self, tree):
        with tree.inputs.panel("Peptide", default_closed=True):
            ala = tree.inputs.color(
                "ALA",
                (0.0, 0.0, 0.0, 1.0),
                description="Set the color for the residue ALA",
            )
            arg = tree.inputs.color(
                "ARG",
                (0.09462588, 0.09462588, 0.09462588, 1.0),
                description="Set the color for the residue ARG",
            )
            asn = tree.inputs.color(
                "ASN",
                (0.11411894, 0.22025453, 0.18504784, 1.0),
                description="Set the color for the residue ASN",
            )
            asp = tree.inputs.color(
                "ASP",
                (0.14087833, 0.3373804, 0.4292262, 1.0),
                description="Set the color for the residue ASP",
            )
            cys = tree.inputs.color(
                "CYS",
                (0.2810098, 0.15101823, 0.10579018, 1.0),
                description="Set the color for the residue CYS",
            )
            glu = tree.inputs.color(
                "GLU",
                (0.20190106, 0.20190106, 0.20190106, 1.0),
                description="Set the color for the residue GLU",
            )
            gln = tree.inputs.color(
                "GLN",
                (0.12770425, 0.20473987, 0.8, 1.0),
                description="Set the color for the residue GLN",
            )
            gly = tree.inputs.color(
                "GLY",
                (0.8, 0.06895567, 0.06778784, 1.0),
                description="Set the color for the residue GLY",
            )
            his = tree.inputs.color(
                "HIS",
                (0.1836438, 0.7612732, 0.341103, 1.0),
                description="Set the color for the residue HIS",
            )
            ile = tree.inputs.color(
                "ILE",
                (0.09119512, 0.6266629, 0.1294339, 1.0),
                description="Set the color for the residue ILE",
            )
            leu = tree.inputs.color(
                "LEU",
                (0.03669303, 0.17026514, 0.41093895, 1.0),
                description="Set the color for the residue LEU",
            )
            lys = tree.inputs.color(
                "LYS",
                (0.05221813, 0.05221813, 0.05221813, 1.0),
                description="Set the color for the residue LYS",
            )
            met = tree.inputs.color(
                "MET",
                (0.5277193, 0.45242295, 0.4978911, 1.0),
                description="Set the color for the residue MET",
            )
            phe = tree.inputs.color(
                "PHE",
                (0.358841, 0.3051302, 0.09418508, 1.0),
                description="Set the color for the residue PHE",
            )
            pro = tree.inputs.color(
                "PRO",
                (0.8, 0.17181273, 0.5252497, 1.0),
                description="Set the color for the residue PRO",
            )
            ser = tree.inputs.color(
                "SER",
                (0.8, 0.722058, 0.05199071, 1.0),
                description="Set the color for the residue SER",
            )
            thr = tree.inputs.color(
                "THR",
                (0.10636823, 1.0, 0.11561122, 1.0),
                description="Set the color for the residue THR",
            )
            trp = tree.inputs.color(
                "TRP",
                (0.5277193, 0.1365669, 0.4106138, 1.0),
                description="Set the color for the residue TRP",
            )
            tyr = tree.inputs.color(
                "TYR",
                (0.08601969, 0.3647821, 0.638272, 1.0),
                description="Set the color for the residue TYR",
            )
            val = tree.inputs.color(
                "VAL",
                (0.0692923, 0.15196387, 0.5596004, 1.0),
                description="Set the color for the residue VAL",
            )
        with tree.inputs.panel("Nucleic", default_closed=True):
            a = tree.inputs.color(
                "A",
                (0.273778, 0.547823, 0.8, 1.0),
                description="Set the color for the residue A",
            )
            c_ = tree.inputs.color(
                "C",
                (0.294582, 0.8, 0.187789, 1.0),
                description="Set the color for the residue C",
            )
            g_ = tree.inputs.color(
                "G",
                (0.85, 0.2514024, 0.17788057, 1.0),
                description="Set the color for the residue G",
            )
            t = tree.inputs.color(
                "T",
                (0.8, 0.269803, 0.526898, 1.0),
                description="Set the color for the residue T",
            )
            ra = tree.inputs.color(
                "rA",
                (0.273778, 0.547823, 0.8, 1.0),
                description="Set the color for the residue rA",
            )
            rc = tree.inputs.color(
                "rC",
                (0.294582, 0.8, 0.187789, 1.0),
                description="Set the color for the residue rC",
            )
            rg = tree.inputs.color(
                "rG",
                (0.85, 0.2514024, 0.17788057, 1.0),
                description="Set the color for the residue rG",
            )
            ru = tree.inputs.color(
                "rU",
                (0.8, 0.269803, 0.526898, 1.0),
                description="Set the color for the residue rU",
            )
        color = tree.outputs.color("Color", (0.8, 0.8, 0.8, 1.0))

        group = ResidueName()
        index_switch = g.IndexSwitch.color(
            group,
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
                (0.8, 0.8, 0.8, 1.0),
                (0.8, 0.8, 0.8, 1.0),
                (0.8, 0.8, 0.8, 1.0),
                (0.8, 0.8, 0.8, 1.0),
                (0.8, 0.8, 0.8, 1.0),
                (0.8, 0.8, 0.8, 1.0),
                (0.8, 0.8, 0.8, 1.0),
                (0.8, 0.8, 0.8, 1.0),
                (0.8, 0.8, 0.8, 1.0),
                (0.8, 0.8, 0.8, 1.0),
                a,
                c_,
                g_,
                t,
                (0.8, 0.8, 0.8, 1.0),
                (0.8, 0.8, 0.8, 1.0),
                (0.8, 0.8, 0.8, 1.0),
                (0.8, 0.8, 0.8, 1.0),
                (0.8, 0.8, 0.8, 1.0),
                (0.8, 0.8, 0.8, 1.0),
                ra,
                rc,
                rg,
                ru,
            ),
        )
        (
            BetweenInteger(value=group, upper=43).o.boolean.switch.color(
                Color(), index_switch
            )
            >> color
        )


ASSET = ColorResName

ASSET_METADATA = {
    "catalog_id": "d3f975df-8408-4972-a669-8187a57e01d0",
}
