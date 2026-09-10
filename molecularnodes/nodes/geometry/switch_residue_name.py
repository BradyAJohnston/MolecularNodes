# Node-group asset "Switch Residue Name" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputInteger
from .residue_name import ResidueName


class SwitchResidueName(AssetGeometryGroup):
    """
    Switch Residue Name

    Parameters
    ----------
    ala : InputInteger
        ALA
    arg : InputInteger
        ARG
    asn : InputInteger
        ASN
    asp : InputInteger
        ASP
    cys : InputInteger
        CYS
    glu : InputInteger
        GLU
    gln : InputInteger
        GLN
    gly : InputInteger
        GLY
    his : InputInteger
        HIS
    ile : InputInteger
        ILE
    leu : InputInteger
        LEU
    lys : InputInteger
        LYS
    met : InputInteger
        MET
    phe : InputInteger
        PHE
    pro : InputInteger
        PRO
    ser : InputInteger
        SER
    thr : InputInteger
        THR
    trp : InputInteger
        TRP
    tyr : InputInteger
        TYR
    val : InputInteger
        VAL

    Inputs
    ------
    i.ala : IntegerSocket
        ALA
    i.arg : IntegerSocket
        ARG
    i.asn : IntegerSocket
        ASN
    i.asp : IntegerSocket
        ASP
    i.cys : IntegerSocket
        CYS
    i.glu : IntegerSocket
        GLU
    i.gln : IntegerSocket
        GLN
    i.gly : IntegerSocket
        GLY
    i.his : IntegerSocket
        HIS
    i.ile : IntegerSocket
        ILE
    i.leu : IntegerSocket
        LEU
    i.lys : IntegerSocket
        LYS
    i.met : IntegerSocket
        MET
    i.phe : IntegerSocket
        PHE
    i.pro : IntegerSocket
        PRO
    i.ser : IntegerSocket
        SER
    i.thr : IntegerSocket
        THR
    i.trp : IntegerSocket
        TRP
    i.tyr : IntegerSocket
        TYR
    i.val : IntegerSocket
        VAL

    Outputs
    -------
    o.output : IntegerSocket
        Output
    """

    _name = "Switch Residue Name"
    _asset_name = "Switch Residue Name"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"

    class _Inputs(SocketAccessor):
        ala: IntegerSocket
        """ALA"""
        arg: IntegerSocket
        """ARG"""
        asn: IntegerSocket
        """ASN"""
        asp: IntegerSocket
        """ASP"""
        cys: IntegerSocket
        """CYS"""
        glu: IntegerSocket
        """GLU"""
        gln: IntegerSocket
        """GLN"""
        gly: IntegerSocket
        """GLY"""
        his: IntegerSocket
        """HIS"""
        ile: IntegerSocket
        """ILE"""
        leu: IntegerSocket
        """LEU"""
        lys: IntegerSocket
        """LYS"""
        met: IntegerSocket
        """MET"""
        phe: IntegerSocket
        """PHE"""
        pro: IntegerSocket
        """PRO"""
        ser: IntegerSocket
        """SER"""
        thr: IntegerSocket
        """THR"""
        trp: IntegerSocket
        """TRP"""
        tyr: IntegerSocket
        """TYR"""
        val: IntegerSocket
        """VAL"""

    class _Outputs(SocketAccessor):
        output: IntegerSocket
        """Output"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        ala: InputInteger = 6,
        arg: InputInteger = 6,
        asn: InputInteger = 6,
        asp: InputInteger = 6,
        cys: InputInteger = 6,
        glu: InputInteger = 6,
        gln: InputInteger = 6,
        gly: InputInteger = 6,
        his: InputInteger = 6,
        ile: InputInteger = 6,
        leu: InputInteger = 6,
        lys: InputInteger = 6,
        met: InputInteger = 6,
        phe: InputInteger = 6,
        pro: InputInteger = 6,
        ser: InputInteger = 6,
        thr: InputInteger = 6,
        trp: InputInteger = 6,
        tyr: InputInteger = 6,
        val: InputInteger = 6,
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

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        ala = tree.inputs.integer("ALA", 6)
        arg = tree.inputs.integer("ARG", 6)
        asn = tree.inputs.integer("ASN", 6)
        asp = tree.inputs.integer("ASP", 6)
        cys = tree.inputs.integer("CYS", 6)
        glu = tree.inputs.integer("GLU", 6)
        gln = tree.inputs.integer("GLN", 6)
        gly = tree.inputs.integer("GLY", 6)
        his = tree.inputs.integer("HIS", 6)
        ile = tree.inputs.integer("ILE", 6)
        leu = tree.inputs.integer("LEU", 6)
        lys = tree.inputs.integer("LYS", 6)
        met = tree.inputs.integer("MET", 6)
        phe = tree.inputs.integer("PHE", 6)
        pro = tree.inputs.integer("PRO", 6)
        ser = tree.inputs.integer("SER", 6)
        thr = tree.inputs.integer("THR", 6)
        trp = tree.inputs.integer("TRP", 6)
        tyr = tree.inputs.integer("TYR", 6)
        val = tree.inputs.integer("VAL", 6)
        output = tree.outputs.integer("Output")

        (
            g.IndexSwitch.integer(
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
            >> output
        )


ASSET = SwitchResidueName

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
