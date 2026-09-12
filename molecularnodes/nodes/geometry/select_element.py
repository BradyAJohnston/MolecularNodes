# Node-group asset "Select Element" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from .atomic_number import AtomicNumber
from .boolean_andor import BooleanAndOr


class SelectElement(AssetGeometryGroup):
    """
    Select Element

    Parameters
    ----------
    and_ : InputBoolean
        The resulting selection must overlap with this input selection
    or_ : InputBoolean
        The resulting selection can be calculated from this node or be from this input selection
    h : InputBoolean
        Select the element H
    he : InputBoolean
        Select the element He
    li : InputBoolean
        Select the element Li
    be : InputBoolean
        Select the element Be
    b : InputBoolean
        Select the element B
    c : InputBoolean
        Select the element C
    n : InputBoolean
        Select the element N
    o : InputBoolean
        Select the element O
    f : InputBoolean
        Select the element F
    ne : InputBoolean
        Select the element Ne
    na : InputBoolean
        Select the element Na
    mg : InputBoolean
        Select the element Mg
    al : InputBoolean
        Select the element Al
    si : InputBoolean
        Select the element Si
    p : InputBoolean
        Select the element P
    s : InputBoolean
        Select the element S
    cl : InputBoolean
        Select the element Cl
    ar : InputBoolean
        Select the element Ar
    k : InputBoolean
        Select the element K
    ca : InputBoolean
        Select the element Ca
    sc : InputBoolean
        Select the element Sc
    ti : InputBoolean
        Select the element Ti
    v : InputBoolean
        Select the element V
    cr : InputBoolean
        Select the element Cr
    mn : InputBoolean
        Select the element Mn
    fe : InputBoolean
        Select the element Fe
    co : InputBoolean
        Select the element Co
    ni : InputBoolean
        Select the element Ni
    cu : InputBoolean
        Select the element Cu
    zn : InputBoolean
        Select the element Zn
    ga : InputBoolean
        Select the element Ga
    ge : InputBoolean
        Select the element Ge
    as_ : InputBoolean
        Select the element As
    se : InputBoolean
        Select the element Se
    br : InputBoolean
        Select the element Br
    kr : InputBoolean
        Select the element Kr
    rb : InputBoolean
        Select the element Rb
    sr : InputBoolean
        Select the element Sr
    y : InputBoolean
        Select the element Y
    zr : InputBoolean
        Select the element Zr
    nb : InputBoolean
        Select the element Nb
    mo : InputBoolean
        Select the element Mo
    tc : InputBoolean
        Select the element Tc
    ru : InputBoolean
        Select the element Ru
    rh : InputBoolean
        Select the element Rh
    pd : InputBoolean
        Select the element Pd
    ag : InputBoolean
        Select the element Ag
    cd : InputBoolean
        Select the element Cd
    in_ : InputBoolean
        Select the element In
    sn : InputBoolean
        Select the element Sn
    sb : InputBoolean
        Select the element Sb
    te : InputBoolean
        Select the element Te
    i : InputBoolean
        Select the element I
    xe : InputBoolean
        Select the element Xe
    cs : InputBoolean
        Select the element Cs
    ba : InputBoolean
        Select the element Ba
    la : InputBoolean
        Select the element La
    ce : InputBoolean
        Select the element Ce
    pr : InputBoolean
        Select the element Pr
    nd : InputBoolean
        Select the element Nd
    pm : InputBoolean
        Select the element Pm
    sm : InputBoolean
        Select the element Sm
    eu : InputBoolean
        Select the element Eu
    gd : InputBoolean
        Select the element Gd
    tb : InputBoolean
        Select the element Tb
    dy : InputBoolean
        Select the element Dy
    ho : InputBoolean
        Select the element Ho
    er : InputBoolean
        Select the element Er
    tm : InputBoolean
        Select the element Tm
    yb : InputBoolean
        Select the element Yb
    lu : InputBoolean
        Select the element Lu
    hf : InputBoolean
        Select the element Hf
    ta : InputBoolean
        Select the element Ta
    w : InputBoolean
        Select the element W
    re : InputBoolean
        Select the element Re
    os : InputBoolean
        Select the element Os
    ir : InputBoolean
        Select the element Ir
    pt : InputBoolean
        Select the element Pt
    au : InputBoolean
        Select the element Au
    hg : InputBoolean
        Select the element Hg

    Inputs
    ------
    i.and_ : BooleanSocket
        The resulting selection must overlap with this input selection
    i.or_ : BooleanSocket
        The resulting selection can be calculated from this node or be from this input selection
    i.h : BooleanSocket
        Select the element H
    i.he : BooleanSocket
        Select the element He
    i.li : BooleanSocket
        Select the element Li
    i.be : BooleanSocket
        Select the element Be
    i.b : BooleanSocket
        Select the element B
    i.c : BooleanSocket
        Select the element C
    i.n : BooleanSocket
        Select the element N
    i.o : BooleanSocket
        Select the element O
    i.f : BooleanSocket
        Select the element F
    i.ne : BooleanSocket
        Select the element Ne
    i.na : BooleanSocket
        Select the element Na
    i.mg : BooleanSocket
        Select the element Mg
    i.al : BooleanSocket
        Select the element Al
    i.si : BooleanSocket
        Select the element Si
    i.p : BooleanSocket
        Select the element P
    i.s : BooleanSocket
        Select the element S
    i.cl : BooleanSocket
        Select the element Cl
    i.ar : BooleanSocket
        Select the element Ar
    i.k : BooleanSocket
        Select the element K
    i.ca : BooleanSocket
        Select the element Ca
    i.sc : BooleanSocket
        Select the element Sc
    i.ti : BooleanSocket
        Select the element Ti
    i.v : BooleanSocket
        Select the element V
    i.cr : BooleanSocket
        Select the element Cr
    i.mn : BooleanSocket
        Select the element Mn
    i.fe : BooleanSocket
        Select the element Fe
    i.co : BooleanSocket
        Select the element Co
    i.ni : BooleanSocket
        Select the element Ni
    i.cu : BooleanSocket
        Select the element Cu
    i.zn : BooleanSocket
        Select the element Zn
    i.ga : BooleanSocket
        Select the element Ga
    i.ge : BooleanSocket
        Select the element Ge
    i.as_ : BooleanSocket
        Select the element As
    i.se : BooleanSocket
        Select the element Se
    i.br : BooleanSocket
        Select the element Br
    i.kr : BooleanSocket
        Select the element Kr
    i.rb : BooleanSocket
        Select the element Rb
    i.sr : BooleanSocket
        Select the element Sr
    i.y : BooleanSocket
        Select the element Y
    i.zr : BooleanSocket
        Select the element Zr
    i.nb : BooleanSocket
        Select the element Nb
    i.mo : BooleanSocket
        Select the element Mo
    i.tc : BooleanSocket
        Select the element Tc
    i.ru : BooleanSocket
        Select the element Ru
    i.rh : BooleanSocket
        Select the element Rh
    i.pd : BooleanSocket
        Select the element Pd
    i.ag : BooleanSocket
        Select the element Ag
    i.cd : BooleanSocket
        Select the element Cd
    i.in_ : BooleanSocket
        Select the element In
    i.sn : BooleanSocket
        Select the element Sn
    i.sb : BooleanSocket
        Select the element Sb
    i.te : BooleanSocket
        Select the element Te
    i.i : BooleanSocket
        Select the element I
    i.xe : BooleanSocket
        Select the element Xe
    i.cs : BooleanSocket
        Select the element Cs
    i.ba : BooleanSocket
        Select the element Ba
    i.la : BooleanSocket
        Select the element La
    i.ce : BooleanSocket
        Select the element Ce
    i.pr : BooleanSocket
        Select the element Pr
    i.nd : BooleanSocket
        Select the element Nd
    i.pm : BooleanSocket
        Select the element Pm
    i.sm : BooleanSocket
        Select the element Sm
    i.eu : BooleanSocket
        Select the element Eu
    i.gd : BooleanSocket
        Select the element Gd
    i.tb : BooleanSocket
        Select the element Tb
    i.dy : BooleanSocket
        Select the element Dy
    i.ho : BooleanSocket
        Select the element Ho
    i.er : BooleanSocket
        Select the element Er
    i.tm : BooleanSocket
        Select the element Tm
    i.yb : BooleanSocket
        Select the element Yb
    i.lu : BooleanSocket
        Select the element Lu
    i.hf : BooleanSocket
        Select the element Hf
    i.ta : BooleanSocket
        Select the element Ta
    i.w : BooleanSocket
        Select the element W
    i.re : BooleanSocket
        Select the element Re
    i.os : BooleanSocket
        Select the element Os
    i.ir : BooleanSocket
        Select the element Ir
    i.pt : BooleanSocket
        Select the element Pt
    i.au : BooleanSocket
        Select the element Au
    i.hg : BooleanSocket
        Select the element Hg

    Outputs
    -------
    o.selection : BooleanSocket
        The calculated selection
    o.inverted : BooleanSocket
        The inverse of the calculated selection
    """

    _name = "Select Element"
    _asset_name = "Select Element"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.select_element"}

    class _Inputs(SocketAccessor):
        and_: BooleanSocket
        """The resulting selection must overlap with this input selection"""
        or_: BooleanSocket
        """The resulting selection can be calculated from this node or be from this input selection"""
        h: BooleanSocket
        """Select the element H"""
        he: BooleanSocket
        """Select the element He"""
        li: BooleanSocket
        """Select the element Li"""
        be: BooleanSocket
        """Select the element Be"""
        b: BooleanSocket
        """Select the element B"""
        c: BooleanSocket
        """Select the element C"""
        n: BooleanSocket
        """Select the element N"""
        o: BooleanSocket
        """Select the element O"""
        f: BooleanSocket
        """Select the element F"""
        ne: BooleanSocket
        """Select the element Ne"""
        na: BooleanSocket
        """Select the element Na"""
        mg: BooleanSocket
        """Select the element Mg"""
        al: BooleanSocket
        """Select the element Al"""
        si: BooleanSocket
        """Select the element Si"""
        p: BooleanSocket
        """Select the element P"""
        s: BooleanSocket
        """Select the element S"""
        cl: BooleanSocket
        """Select the element Cl"""
        ar: BooleanSocket
        """Select the element Ar"""
        k: BooleanSocket
        """Select the element K"""
        ca: BooleanSocket
        """Select the element Ca"""
        sc: BooleanSocket
        """Select the element Sc"""
        ti: BooleanSocket
        """Select the element Ti"""
        v: BooleanSocket
        """Select the element V"""
        cr: BooleanSocket
        """Select the element Cr"""
        mn: BooleanSocket
        """Select the element Mn"""
        fe: BooleanSocket
        """Select the element Fe"""
        co: BooleanSocket
        """Select the element Co"""
        ni: BooleanSocket
        """Select the element Ni"""
        cu: BooleanSocket
        """Select the element Cu"""
        zn: BooleanSocket
        """Select the element Zn"""
        ga: BooleanSocket
        """Select the element Ga"""
        ge: BooleanSocket
        """Select the element Ge"""
        as_: BooleanSocket
        """Select the element As"""
        se: BooleanSocket
        """Select the element Se"""
        br: BooleanSocket
        """Select the element Br"""
        kr: BooleanSocket
        """Select the element Kr"""
        rb: BooleanSocket
        """Select the element Rb"""
        sr: BooleanSocket
        """Select the element Sr"""
        y: BooleanSocket
        """Select the element Y"""
        zr: BooleanSocket
        """Select the element Zr"""
        nb: BooleanSocket
        """Select the element Nb"""
        mo: BooleanSocket
        """Select the element Mo"""
        tc: BooleanSocket
        """Select the element Tc"""
        ru: BooleanSocket
        """Select the element Ru"""
        rh: BooleanSocket
        """Select the element Rh"""
        pd: BooleanSocket
        """Select the element Pd"""
        ag: BooleanSocket
        """Select the element Ag"""
        cd: BooleanSocket
        """Select the element Cd"""
        in_: BooleanSocket
        """Select the element In"""
        sn: BooleanSocket
        """Select the element Sn"""
        sb: BooleanSocket
        """Select the element Sb"""
        te: BooleanSocket
        """Select the element Te"""
        i: BooleanSocket
        """Select the element I"""
        xe: BooleanSocket
        """Select the element Xe"""
        cs: BooleanSocket
        """Select the element Cs"""
        ba: BooleanSocket
        """Select the element Ba"""
        la: BooleanSocket
        """Select the element La"""
        ce: BooleanSocket
        """Select the element Ce"""
        pr: BooleanSocket
        """Select the element Pr"""
        nd: BooleanSocket
        """Select the element Nd"""
        pm: BooleanSocket
        """Select the element Pm"""
        sm: BooleanSocket
        """Select the element Sm"""
        eu: BooleanSocket
        """Select the element Eu"""
        gd: BooleanSocket
        """Select the element Gd"""
        tb: BooleanSocket
        """Select the element Tb"""
        dy: BooleanSocket
        """Select the element Dy"""
        ho: BooleanSocket
        """Select the element Ho"""
        er: BooleanSocket
        """Select the element Er"""
        tm: BooleanSocket
        """Select the element Tm"""
        yb: BooleanSocket
        """Select the element Yb"""
        lu: BooleanSocket
        """Select the element Lu"""
        hf: BooleanSocket
        """Select the element Hf"""
        ta: BooleanSocket
        """Select the element Ta"""
        w: BooleanSocket
        """Select the element W"""
        re: BooleanSocket
        """Select the element Re"""
        os: BooleanSocket
        """Select the element Os"""
        ir: BooleanSocket
        """Select the element Ir"""
        pt: BooleanSocket
        """Select the element Pt"""
        au: BooleanSocket
        """Select the element Au"""
        hg: BooleanSocket
        """Select the element Hg"""

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
        h: InputBoolean = False,
        he: InputBoolean = False,
        li: InputBoolean = False,
        be: InputBoolean = False,
        b: InputBoolean = False,
        c: InputBoolean = False,
        n: InputBoolean = False,
        o: InputBoolean = False,
        f: InputBoolean = False,
        ne: InputBoolean = False,
        na: InputBoolean = False,
        mg: InputBoolean = False,
        al: InputBoolean = False,
        si: InputBoolean = False,
        p: InputBoolean = False,
        s: InputBoolean = False,
        cl: InputBoolean = False,
        ar: InputBoolean = False,
        k: InputBoolean = False,
        ca: InputBoolean = False,
        sc: InputBoolean = False,
        ti: InputBoolean = False,
        v: InputBoolean = False,
        cr: InputBoolean = False,
        mn: InputBoolean = False,
        fe: InputBoolean = False,
        co: InputBoolean = False,
        ni: InputBoolean = False,
        cu: InputBoolean = False,
        zn: InputBoolean = False,
        ga: InputBoolean = False,
        ge: InputBoolean = False,
        as_: InputBoolean = False,
        se: InputBoolean = False,
        br: InputBoolean = False,
        kr: InputBoolean = False,
        rb: InputBoolean = False,
        sr: InputBoolean = False,
        y: InputBoolean = False,
        zr: InputBoolean = False,
        nb: InputBoolean = False,
        mo: InputBoolean = False,
        tc: InputBoolean = False,
        ru: InputBoolean = False,
        rh: InputBoolean = False,
        pd: InputBoolean = False,
        ag: InputBoolean = False,
        cd: InputBoolean = False,
        in_: InputBoolean = False,
        sn: InputBoolean = False,
        sb: InputBoolean = False,
        te: InputBoolean = False,
        i: InputBoolean = False,
        xe: InputBoolean = False,
        cs: InputBoolean = False,
        ba: InputBoolean = False,
        la: InputBoolean = False,
        ce: InputBoolean = False,
        pr: InputBoolean = False,
        nd: InputBoolean = False,
        pm: InputBoolean = False,
        sm: InputBoolean = False,
        eu: InputBoolean = False,
        gd: InputBoolean = False,
        tb: InputBoolean = False,
        dy: InputBoolean = False,
        ho: InputBoolean = False,
        er: InputBoolean = False,
        tm: InputBoolean = False,
        yb: InputBoolean = False,
        lu: InputBoolean = False,
        hf: InputBoolean = False,
        ta: InputBoolean = False,
        w: InputBoolean = False,
        re: InputBoolean = False,
        os: InputBoolean = False,
        ir: InputBoolean = False,
        pt: InputBoolean = False,
        au: InputBoolean = False,
        hg: InputBoolean = False,
    ):
        super().__init__(
            **{
                "And": and_,
                "Or": or_,
                "H": h,
                "He": he,
                "Li": li,
                "Be": be,
                "B": b,
                "C": c,
                "N": n,
                "O": o,
                "F": f,
                "Ne": ne,
                "Na": na,
                "Mg": mg,
                "Al": al,
                "Si": si,
                "P": p,
                "S": s,
                "Cl": cl,
                "Ar": ar,
                "K": k,
                "Ca": ca,
                "Sc": sc,
                "Ti": ti,
                "V": v,
                "Cr": cr,
                "Mn": mn,
                "Fe": fe,
                "Co": co,
                "Ni": ni,
                "Cu": cu,
                "Zn": zn,
                "Ga": ga,
                "Ge": ge,
                "As": as_,
                "Se": se,
                "Br": br,
                "Kr": kr,
                "Rb": rb,
                "Sr": sr,
                "Y": y,
                "Zr": zr,
                "Nb": nb,
                "Mo": mo,
                "Tc": tc,
                "Ru": ru,
                "Rh": rh,
                "Pd": pd,
                "Ag": ag,
                "Cd": cd,
                "In": in_,
                "Sn": sn,
                "Sb": sb,
                "Te": te,
                "I": i,
                "Xe": xe,
                "Cs": cs,
                "Ba": ba,
                "La": la,
                "Ce": ce,
                "Pr": pr,
                "Nd": nd,
                "Pm": pm,
                "Sm": sm,
                "Eu": eu,
                "Gd": gd,
                "Tb": tb,
                "Dy": dy,
                "Ho": ho,
                "Er": er,
                "Tm": tm,
                "Yb": yb,
                "Lu": lu,
                "Hf": hf,
                "Ta": ta,
                "W": w,
                "Re": re,
                "Os": os,
                "Ir": ir,
                "Pt": pt,
                "Au": au,
                "Hg": hg,
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
        with tree.inputs.panel("1-20", default_closed=True):
            h = tree.inputs.boolean("H", False, description="Select the element H")
            he = tree.inputs.boolean("He", False, description="Select the element He")
            li = tree.inputs.boolean("Li", False, description="Select the element Li")
            be = tree.inputs.boolean("Be", False, description="Select the element Be")
            b = tree.inputs.boolean("B", False, description="Select the element B")
            c_ = tree.inputs.boolean("C", False, description="Select the element C")
            n = tree.inputs.boolean("N", False, description="Select the element N")
            o = tree.inputs.boolean("O", False, description="Select the element O")
            f = tree.inputs.boolean("F", False, description="Select the element F")
            ne = tree.inputs.boolean("Ne", False, description="Select the element Ne")
            na = tree.inputs.boolean("Na", False, description="Select the element Na")
            mg = tree.inputs.boolean("Mg", False, description="Select the element Mg")
            al = tree.inputs.boolean("Al", False, description="Select the element Al")
            si = tree.inputs.boolean("Si", False, description="Select the element Si")
            p = tree.inputs.boolean("P", False, description="Select the element P")
            s_ = tree.inputs.boolean("S", False, description="Select the element S")
            cl = tree.inputs.boolean("Cl", False, description="Select the element Cl")
            ar = tree.inputs.boolean("Ar", False, description="Select the element Ar")
            k = tree.inputs.boolean("K", False, description="Select the element K")
            ca = tree.inputs.boolean("Ca", False, description="Select the element Ca")
        with tree.inputs.panel("21-40", default_closed=True):
            sc = tree.inputs.boolean("Sc", False, description="Select the element Sc")
            ti = tree.inputs.boolean("Ti", False, description="Select the element Ti")
            v = tree.inputs.boolean("V", False, description="Select the element V")
            cr = tree.inputs.boolean("Cr", False, description="Select the element Cr")
            mn = tree.inputs.boolean("Mn", False, description="Select the element Mn")
            fe = tree.inputs.boolean("Fe", False, description="Select the element Fe")
            co = tree.inputs.boolean("Co", False, description="Select the element Co")
            ni = tree.inputs.boolean("Ni", False, description="Select the element Ni")
            cu = tree.inputs.boolean("Cu", False, description="Select the element Cu")
            zn = tree.inputs.boolean("Zn", False, description="Select the element Zn")
            ga = tree.inputs.boolean("Ga", False, description="Select the element Ga")
            ge = tree.inputs.boolean("Ge", False, description="Select the element Ge")
            as_ = tree.inputs.boolean("As", False, description="Select the element As")
            se = tree.inputs.boolean("Se", False, description="Select the element Se")
            br = tree.inputs.boolean("Br", False, description="Select the element Br")
            kr = tree.inputs.boolean("Kr", False, description="Select the element Kr")
            rb = tree.inputs.boolean("Rb", False, description="Select the element Rb")
            sr = tree.inputs.boolean("Sr", False, description="Select the element Sr")
            y = tree.inputs.boolean("Y", False, description="Select the element Y")
            zr = tree.inputs.boolean("Zr", False, description="Select the element Zr")
        with tree.inputs.panel("41-60", default_closed=True):
            nb = tree.inputs.boolean("Nb", False, description="Select the element Nb")
            mo = tree.inputs.boolean("Mo", False, description="Select the element Mo")
            tc = tree.inputs.boolean("Tc", False, description="Select the element Tc")
            ru = tree.inputs.boolean("Ru", False, description="Select the element Ru")
            rh = tree.inputs.boolean("Rh", False, description="Select the element Rh")
            pd = tree.inputs.boolean("Pd", False, description="Select the element Pd")
            ag = tree.inputs.boolean("Ag", False, description="Select the element Ag")
            cd = tree.inputs.boolean("Cd", False, description="Select the element Cd")
            in_ = tree.inputs.boolean("In", False, description="Select the element In")
            sn = tree.inputs.boolean("Sn", False, description="Select the element Sn")
            sb = tree.inputs.boolean("Sb", False, description="Select the element Sb")
            te = tree.inputs.boolean("Te", False, description="Select the element Te")
            i = tree.inputs.boolean("I", False, description="Select the element I")
            xe = tree.inputs.boolean("Xe", False, description="Select the element Xe")
            cs = tree.inputs.boolean("Cs", False, description="Select the element Cs")
            ba = tree.inputs.boolean("Ba", False, description="Select the element Ba")
            la = tree.inputs.boolean("La", False, description="Select the element La")
            ce = tree.inputs.boolean("Ce", False, description="Select the element Ce")
            pr = tree.inputs.boolean("Pr", False, description="Select the element Pr")
            nd = tree.inputs.boolean("Nd", False, description="Select the element Nd")
        with tree.inputs.panel("61-80", default_closed=True):
            pm = tree.inputs.boolean("Pm", False, description="Select the element Pm")
            sm = tree.inputs.boolean("Sm", False, description="Select the element Sm")
            eu = tree.inputs.boolean("Eu", False, description="Select the element Eu")
            gd = tree.inputs.boolean("Gd", False, description="Select the element Gd")
            tb = tree.inputs.boolean("Tb", False, description="Select the element Tb")
            dy = tree.inputs.boolean("Dy", False, description="Select the element Dy")
            ho = tree.inputs.boolean("Ho", False, description="Select the element Ho")
            er = tree.inputs.boolean("Er", False, description="Select the element Er")
            tm = tree.inputs.boolean("Tm", False, description="Select the element Tm")
            yb = tree.inputs.boolean("Yb", False, description="Select the element Yb")
            lu = tree.inputs.boolean("Lu", False, description="Select the element Lu")
            hf = tree.inputs.boolean("Hf", False, description="Select the element Hf")
            ta = tree.inputs.boolean("Ta", False, description="Select the element Ta")
            w = tree.inputs.boolean("W", False, description="Select the element W")
            re = tree.inputs.boolean("Re", False, description="Select the element Re")
            os = tree.inputs.boolean("Os", False, description="Select the element Os")
            ir = tree.inputs.boolean("Ir", False, description="Select the element Ir")
            pt = tree.inputs.boolean("Pt", False, description="Select the element Pt")
            au = tree.inputs.boolean("Au", False, description="Select the element Au")
            hg = tree.inputs.boolean("Hg", False, description="Select the element Hg")
        selection = tree.outputs.boolean(
            "Selection", description="The calculated selection"
        )
        inverted = tree.outputs.boolean(
            "Inverted", description="The inverse of the calculated selection"
        )

        index_switch = g.IndexSwitch.boolean(
            AtomicNumber(),
            (
                False,
                h,
                he,
                li,
                be,
                b,
                c_,
                n,
                o,
                f,
                ne,
                na,
                mg,
                al,
                si,
                p,
                s_,
                cl,
                ar,
                k,
                ca,
                sc,
                ti,
                v,
                cr,
                mn,
                fe,
                co,
                ni,
                cu,
                zn,
                ga,
                ge,
                as_,
                se,
                br,
                kr,
                rb,
                sr,
                y,
                zr,
                nb,
                mo,
                tc,
                ru,
                rh,
                pd,
                ag,
                cd,
                in_,
                sn,
                sb,
                te,
                i,
                xe,
                cs,
                ba,
                la,
                ce,
                pr,
                nd,
                pm,
                sm,
                eu,
                gd,
                tb,
                dy,
                ho,
                er,
                tm,
                yb,
                lu,
                hf,
                ta,
                w,
                re,
                os,
                ir,
                pt,
                au,
                hg,
            ),
        )
        group = BooleanAndOr(and_=and_, or_=or_, boolean=index_switch)

        group >> selection
        group.o.inverted >> inverted


ASSET = SelectElement

ASSET_METADATA = {
    "catalog_id": "bd1f205b-fea5-4700-b2c2-754f3321e969",
}
