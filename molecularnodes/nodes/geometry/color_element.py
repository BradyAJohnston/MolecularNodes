# Node-group asset "Color Element" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    ColorSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputColor
from .atomic_number import AtomicNumber


class ColorElement(AssetGeometryGroup):
    """
    Color Element

    Parameters
    ----------
    h : InputColor
        Set the color for the element H
    he : InputColor
        Set the color for the element He
    li : InputColor
        Set the color for the element Li
    be : InputColor
        Set the color for the element Be
    b : InputColor
        Set the color for the element B
    c : InputColor
        Set the color for the element C
    n : InputColor
        Set the color for the element N
    o : InputColor
        Set the color for the element O
    f : InputColor
        Set the color for the element F
    ne : InputColor
        Set the color for the element Ne
    na : InputColor
        Set the color for the element Na
    mg : InputColor
        Set the color for the element Mg
    al : InputColor
        Set the color for the element Al
    si : InputColor
        Set the color for the element Si
    p : InputColor
        Set the color for the element P
    s : InputColor
        Set the color for the element S
    cl : InputColor
        Set the color for the element Cl
    ar : InputColor
        Set the color for the element Ar
    k : InputColor
        Set the color for the element K
    ca : InputColor
        Set the color for the element Ca
    sc : InputColor
        Set the color for the element Sc
    ti : InputColor
        Set the color for the element Ti
    v : InputColor
        Set the color for the element V
    cr : InputColor
        Set the color for the element Cr
    mn : InputColor
        Set the color for the element Mn
    fe : InputColor
        Set the color for the element Fe
    co : InputColor
        Set the color for the element Co
    ni : InputColor
        Set the color for the element Ni
    cu : InputColor
        Set the color for the element Cu
    zn : InputColor
        Set the color for the element Zn
    ga : InputColor
        Set the color for the element Ga
    ge : InputColor
        Set the color for the element Ge
    as_ : InputColor
        Set the color for the element As
    se : InputColor
        Set the color for the element Se
    br : InputColor
        Set the color for the element Br
    kr : InputColor
        Set the color for the element Kr
    rb : InputColor
        Set the color for the element Rb
    sr : InputColor
        Set the color for the element Sr
    y : InputColor
        Set the color for the element Y
    zr : InputColor
        Set the color for the element Zr
    nb : InputColor
        Set the color for the element Nb
    mo : InputColor
        Set the color for the element Mo
    tc : InputColor
        Set the color for the element Tc
    ru : InputColor
        Set the color for the element Ru
    rh : InputColor
        Set the color for the element Rh
    pd : InputColor
        Set the color for the element Pd
    ag : InputColor
        Set the color for the element Ag
    cd : InputColor
        Set the color for the element Cd
    in_ : InputColor
        Set the color for the element In
    sn : InputColor
        Set the color for the element Sn
    sb : InputColor
        Set the color for the element Sb
    te : InputColor
        Set the color for the element Te
    i : InputColor
        Set the color for the element I
    xe : InputColor
        Set the color for the element Xe
    cs : InputColor
        Set the color for the element Cs
    ba : InputColor
        Set the color for the element Ba
    la : InputColor
        Set the color for the element La
    ce : InputColor
        Set the color for the element Ce
    pr : InputColor
        Set the color for the element Pr
    nd : InputColor
        Set the color for the element Nd
    pm : InputColor
        Set the color for the element Pm
    sm : InputColor
        Set the color for the element Sm
    eu : InputColor
        Set the color for the element Eu
    gd : InputColor
        Set the color for the element Gd
    tb : InputColor
        Set the color for the element Tb
    dy : InputColor
        Set the color for the element Dy
    ho : InputColor
        Set the color for the element Ho
    er : InputColor
        Set the color for the element Er
    tm : InputColor
        Set the color for the element Tm
    yb : InputColor
        Set the color for the element Yb
    lu : InputColor
        Set the color for the element Lu
    hf : InputColor
        Set the color for the element Hf
    ta : InputColor
        Set the color for the element Ta
    w : InputColor
        Set the color for the element W
    re : InputColor
        Set the color for the element Re
    os : InputColor
        Set the color for the element Os
    ir : InputColor
        Set the color for the element Ir
    pt : InputColor
        Set the color for the element Pt
    au : InputColor
        Set the color for the element Au
    hg : InputColor
        Set the color for the element Hg

    Inputs
    ------
    i.h : ColorSocket
        Set the color for the element H
    i.he : ColorSocket
        Set the color for the element He
    i.li : ColorSocket
        Set the color for the element Li
    i.be : ColorSocket
        Set the color for the element Be
    i.b : ColorSocket
        Set the color for the element B
    i.c : ColorSocket
        Set the color for the element C
    i.n : ColorSocket
        Set the color for the element N
    i.o : ColorSocket
        Set the color for the element O
    i.f : ColorSocket
        Set the color for the element F
    i.ne : ColorSocket
        Set the color for the element Ne
    i.na : ColorSocket
        Set the color for the element Na
    i.mg : ColorSocket
        Set the color for the element Mg
    i.al : ColorSocket
        Set the color for the element Al
    i.si : ColorSocket
        Set the color for the element Si
    i.p : ColorSocket
        Set the color for the element P
    i.s : ColorSocket
        Set the color for the element S
    i.cl : ColorSocket
        Set the color for the element Cl
    i.ar : ColorSocket
        Set the color for the element Ar
    i.k : ColorSocket
        Set the color for the element K
    i.ca : ColorSocket
        Set the color for the element Ca
    i.sc : ColorSocket
        Set the color for the element Sc
    i.ti : ColorSocket
        Set the color for the element Ti
    i.v : ColorSocket
        Set the color for the element V
    i.cr : ColorSocket
        Set the color for the element Cr
    i.mn : ColorSocket
        Set the color for the element Mn
    i.fe : ColorSocket
        Set the color for the element Fe
    i.co : ColorSocket
        Set the color for the element Co
    i.ni : ColorSocket
        Set the color for the element Ni
    i.cu : ColorSocket
        Set the color for the element Cu
    i.zn : ColorSocket
        Set the color for the element Zn
    i.ga : ColorSocket
        Set the color for the element Ga
    i.ge : ColorSocket
        Set the color for the element Ge
    i.as_ : ColorSocket
        Set the color for the element As
    i.se : ColorSocket
        Set the color for the element Se
    i.br : ColorSocket
        Set the color for the element Br
    i.kr : ColorSocket
        Set the color for the element Kr
    i.rb : ColorSocket
        Set the color for the element Rb
    i.sr : ColorSocket
        Set the color for the element Sr
    i.y : ColorSocket
        Set the color for the element Y
    i.zr : ColorSocket
        Set the color for the element Zr
    i.nb : ColorSocket
        Set the color for the element Nb
    i.mo : ColorSocket
        Set the color for the element Mo
    i.tc : ColorSocket
        Set the color for the element Tc
    i.ru : ColorSocket
        Set the color for the element Ru
    i.rh : ColorSocket
        Set the color for the element Rh
    i.pd : ColorSocket
        Set the color for the element Pd
    i.ag : ColorSocket
        Set the color for the element Ag
    i.cd : ColorSocket
        Set the color for the element Cd
    i.in_ : ColorSocket
        Set the color for the element In
    i.sn : ColorSocket
        Set the color for the element Sn
    i.sb : ColorSocket
        Set the color for the element Sb
    i.te : ColorSocket
        Set the color for the element Te
    i.i : ColorSocket
        Set the color for the element I
    i.xe : ColorSocket
        Set the color for the element Xe
    i.cs : ColorSocket
        Set the color for the element Cs
    i.ba : ColorSocket
        Set the color for the element Ba
    i.la : ColorSocket
        Set the color for the element La
    i.ce : ColorSocket
        Set the color for the element Ce
    i.pr : ColorSocket
        Set the color for the element Pr
    i.nd : ColorSocket
        Set the color for the element Nd
    i.pm : ColorSocket
        Set the color for the element Pm
    i.sm : ColorSocket
        Set the color for the element Sm
    i.eu : ColorSocket
        Set the color for the element Eu
    i.gd : ColorSocket
        Set the color for the element Gd
    i.tb : ColorSocket
        Set the color for the element Tb
    i.dy : ColorSocket
        Set the color for the element Dy
    i.ho : ColorSocket
        Set the color for the element Ho
    i.er : ColorSocket
        Set the color for the element Er
    i.tm : ColorSocket
        Set the color for the element Tm
    i.yb : ColorSocket
        Set the color for the element Yb
    i.lu : ColorSocket
        Set the color for the element Lu
    i.hf : ColorSocket
        Set the color for the element Hf
    i.ta : ColorSocket
        Set the color for the element Ta
    i.w : ColorSocket
        Set the color for the element W
    i.re : ColorSocket
        Set the color for the element Re
    i.os : ColorSocket
        Set the color for the element Os
    i.ir : ColorSocket
        Set the color for the element Ir
    i.pt : ColorSocket
        Set the color for the element Pt
    i.au : ColorSocket
        Set the color for the element Au
    i.hg : ColorSocket
        Set the color for the element Hg

    Outputs
    -------
    o.color : ColorSocket
        The selected colors based on the `atomic_number`
    """

    _name = "Color Element"
    _asset_name = "Color Element"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "COLOR"
    _tree_properties = {"node_tool_idname": "geometry.color_element"}

    class _Inputs(SocketAccessor):
        h: ColorSocket
        """Set the color for the element H"""
        he: ColorSocket
        """Set the color for the element He"""
        li: ColorSocket
        """Set the color for the element Li"""
        be: ColorSocket
        """Set the color for the element Be"""
        b: ColorSocket
        """Set the color for the element B"""
        c: ColorSocket
        """Set the color for the element C"""
        n: ColorSocket
        """Set the color for the element N"""
        o: ColorSocket
        """Set the color for the element O"""
        f: ColorSocket
        """Set the color for the element F"""
        ne: ColorSocket
        """Set the color for the element Ne"""
        na: ColorSocket
        """Set the color for the element Na"""
        mg: ColorSocket
        """Set the color for the element Mg"""
        al: ColorSocket
        """Set the color for the element Al"""
        si: ColorSocket
        """Set the color for the element Si"""
        p: ColorSocket
        """Set the color for the element P"""
        s: ColorSocket
        """Set the color for the element S"""
        cl: ColorSocket
        """Set the color for the element Cl"""
        ar: ColorSocket
        """Set the color for the element Ar"""
        k: ColorSocket
        """Set the color for the element K"""
        ca: ColorSocket
        """Set the color for the element Ca"""
        sc: ColorSocket
        """Set the color for the element Sc"""
        ti: ColorSocket
        """Set the color for the element Ti"""
        v: ColorSocket
        """Set the color for the element V"""
        cr: ColorSocket
        """Set the color for the element Cr"""
        mn: ColorSocket
        """Set the color for the element Mn"""
        fe: ColorSocket
        """Set the color for the element Fe"""
        co: ColorSocket
        """Set the color for the element Co"""
        ni: ColorSocket
        """Set the color for the element Ni"""
        cu: ColorSocket
        """Set the color for the element Cu"""
        zn: ColorSocket
        """Set the color for the element Zn"""
        ga: ColorSocket
        """Set the color for the element Ga"""
        ge: ColorSocket
        """Set the color for the element Ge"""
        as_: ColorSocket
        """Set the color for the element As"""
        se: ColorSocket
        """Set the color for the element Se"""
        br: ColorSocket
        """Set the color for the element Br"""
        kr: ColorSocket
        """Set the color for the element Kr"""
        rb: ColorSocket
        """Set the color for the element Rb"""
        sr: ColorSocket
        """Set the color for the element Sr"""
        y: ColorSocket
        """Set the color for the element Y"""
        zr: ColorSocket
        """Set the color for the element Zr"""
        nb: ColorSocket
        """Set the color for the element Nb"""
        mo: ColorSocket
        """Set the color for the element Mo"""
        tc: ColorSocket
        """Set the color for the element Tc"""
        ru: ColorSocket
        """Set the color for the element Ru"""
        rh: ColorSocket
        """Set the color for the element Rh"""
        pd: ColorSocket
        """Set the color for the element Pd"""
        ag: ColorSocket
        """Set the color for the element Ag"""
        cd: ColorSocket
        """Set the color for the element Cd"""
        in_: ColorSocket
        """Set the color for the element In"""
        sn: ColorSocket
        """Set the color for the element Sn"""
        sb: ColorSocket
        """Set the color for the element Sb"""
        te: ColorSocket
        """Set the color for the element Te"""
        i: ColorSocket
        """Set the color for the element I"""
        xe: ColorSocket
        """Set the color for the element Xe"""
        cs: ColorSocket
        """Set the color for the element Cs"""
        ba: ColorSocket
        """Set the color for the element Ba"""
        la: ColorSocket
        """Set the color for the element La"""
        ce: ColorSocket
        """Set the color for the element Ce"""
        pr: ColorSocket
        """Set the color for the element Pr"""
        nd: ColorSocket
        """Set the color for the element Nd"""
        pm: ColorSocket
        """Set the color for the element Pm"""
        sm: ColorSocket
        """Set the color for the element Sm"""
        eu: ColorSocket
        """Set the color for the element Eu"""
        gd: ColorSocket
        """Set the color for the element Gd"""
        tb: ColorSocket
        """Set the color for the element Tb"""
        dy: ColorSocket
        """Set the color for the element Dy"""
        ho: ColorSocket
        """Set the color for the element Ho"""
        er: ColorSocket
        """Set the color for the element Er"""
        tm: ColorSocket
        """Set the color for the element Tm"""
        yb: ColorSocket
        """Set the color for the element Yb"""
        lu: ColorSocket
        """Set the color for the element Lu"""
        hf: ColorSocket
        """Set the color for the element Hf"""
        ta: ColorSocket
        """Set the color for the element Ta"""
        w: ColorSocket
        """Set the color for the element W"""
        re: ColorSocket
        """Set the color for the element Re"""
        os: ColorSocket
        """Set the color for the element Os"""
        ir: ColorSocket
        """Set the color for the element Ir"""
        pt: ColorSocket
        """Set the color for the element Pt"""
        au: ColorSocket
        """Set the color for the element Au"""
        hg: ColorSocket
        """Set the color for the element Hg"""

    class _Outputs(SocketAccessor):
        color: ColorSocket
        """The selected colors based on the `atomic_number`"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        h: InputColor = None,
        he: InputColor = None,
        li: InputColor = None,
        be: InputColor = None,
        b: InputColor = None,
        c: InputColor = None,
        n: InputColor = None,
        o: InputColor = None,
        f: InputColor = None,
        ne: InputColor = None,
        na: InputColor = None,
        mg: InputColor = None,
        al: InputColor = None,
        si: InputColor = None,
        p: InputColor = None,
        s: InputColor = None,
        cl: InputColor = None,
        ar: InputColor = None,
        k: InputColor = None,
        ca: InputColor = None,
        sc: InputColor = None,
        ti: InputColor = None,
        v: InputColor = None,
        cr: InputColor = None,
        mn: InputColor = None,
        fe: InputColor = None,
        co: InputColor = None,
        ni: InputColor = None,
        cu: InputColor = None,
        zn: InputColor = None,
        ga: InputColor = None,
        ge: InputColor = None,
        as_: InputColor = None,
        se: InputColor = None,
        br: InputColor = None,
        kr: InputColor = None,
        rb: InputColor = None,
        sr: InputColor = None,
        y: InputColor = None,
        zr: InputColor = None,
        nb: InputColor = None,
        mo: InputColor = None,
        tc: InputColor = None,
        ru: InputColor = None,
        rh: InputColor = None,
        pd: InputColor = None,
        ag: InputColor = None,
        cd: InputColor = None,
        in_: InputColor = None,
        sn: InputColor = None,
        sb: InputColor = None,
        te: InputColor = None,
        i: InputColor = None,
        xe: InputColor = None,
        cs: InputColor = None,
        ba: InputColor = None,
        la: InputColor = None,
        ce: InputColor = None,
        pr: InputColor = None,
        nd: InputColor = None,
        pm: InputColor = None,
        sm: InputColor = None,
        eu: InputColor = None,
        gd: InputColor = None,
        tb: InputColor = None,
        dy: InputColor = None,
        ho: InputColor = None,
        er: InputColor = None,
        tm: InputColor = None,
        yb: InputColor = None,
        lu: InputColor = None,
        hf: InputColor = None,
        ta: InputColor = None,
        w: InputColor = None,
        re: InputColor = None,
        os: InputColor = None,
        ir: InputColor = None,
        pt: InputColor = None,
        au: InputColor = None,
        hg: InputColor = None,
    ):
        super().__init__(
            **{
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
        with tree.inputs.panel("1-20", default_closed=True):
            h = tree.inputs.color(
                "H", (1.0, 1.0, 1.0, 1.0), description="Set the color for the element H"
            )
            he = tree.inputs.color(
                "He",
                (0.8509804, 1.0, 1.0, 1.0),
                description="Set the color for the element He",
            )
            li = tree.inputs.color(
                "Li",
                (0.8, 0.5019608, 1.0, 1.0),
                description="Set the color for the element Li",
            )
            be = tree.inputs.color(
                "Be",
                (0.7607843, 1.0, 0.0, 1.0),
                description="Set the color for the element Be",
            )
            b = tree.inputs.color(
                "B",
                (1.0, 0.709804, 0.709804, 1.0),
                description="Set the color for the element B",
            )
            c_ = tree.inputs.color(
                "C",
                (0.5647059, 0.5647059, 0.5647059, 1.0),
                description="Set the color for the element C",
            )
            n = tree.inputs.color(
                "N",
                (0.1882353, 0.3137255, 0.972549, 1.0),
                description="Set the color for the element N",
            )
            o = tree.inputs.color(
                "O",
                (1.0, 0.05098039, 0.05098039, 1.0),
                description="Set the color for the element O",
            )
            f = tree.inputs.color(
                "F",
                (0.5647059, 0.8784314, 0.3137255, 1.0),
                description="Set the color for the element F",
            )
            ne = tree.inputs.color(
                "Ne",
                (0.7019608, 0.890196, 0.9607843, 1.0),
                description="Set the color for the element Ne",
            )
            na = tree.inputs.color(
                "Na",
                (0.6705883, 0.3607843, 0.9490196, 1.0),
                description="Set the color for the element Na",
            )
            mg = tree.inputs.color(
                "Mg",
                (0.5411765, 1.0, 0.0, 1.0),
                description="Set the color for the element Mg",
            )
            al = tree.inputs.color(
                "Al",
                (0.7490196, 0.6509804, 0.6509804, 1.0),
                description="Set the color for the element Al",
            )
            si = tree.inputs.color(
                "Si",
                (0.9411765, 0.7843137, 0.627451, 1.0),
                description="Set the color for the element Si",
            )
            p = tree.inputs.color(
                "P",
                (1.0, 0.5019608, 0.0, 1.0),
                description="Set the color for the element P",
            )
            s_ = tree.inputs.color(
                "S",
                (1.0, 1.0, 0.1882353, 1.0),
                description="Set the color for the element S",
            )
            cl = tree.inputs.color(
                "Cl",
                (0.12156863, 0.9411765, 0.12156863, 1.0),
                description="Set the color for the element Cl",
            )
            ar = tree.inputs.color(
                "Ar",
                (0.5019608, 0.8196079, 0.890196, 1.0),
                description="Set the color for the element Ar",
            )
            k = tree.inputs.color(
                "K",
                (0.5607843, 0.2509804, 0.8313726, 1.0),
                description="Set the color for the element K",
            )
            ca = tree.inputs.color(
                "Ca",
                (0.2392157, 1.0, 0.0, 1.0),
                description="Set the color for the element Ca",
            )
        with tree.inputs.panel("21-40", default_closed=True):
            sc = tree.inputs.color(
                "Sc",
                (0.9019608, 0.9019608, 0.9019608, 1.0),
                description="Set the color for the element Sc",
            )
            ti = tree.inputs.color(
                "Ti",
                (0.7490196, 0.7607843, 0.7803922, 1.0),
                description="Set the color for the element Ti",
            )
            v = tree.inputs.color(
                "V",
                (0.6509804, 0.6509804, 0.6705883, 1.0),
                description="Set the color for the element V",
            )
            cr = tree.inputs.color(
                "Cr",
                (0.5411765, 0.6, 0.7803922, 1.0),
                description="Set the color for the element Cr",
            )
            mn = tree.inputs.color(
                "Mn",
                (0.6117647, 0.4784314, 0.7803922, 1.0),
                description="Set the color for the element Mn",
            )
            fe = tree.inputs.color(
                "Fe",
                (0.8784314, 0.4, 0.2, 1.0),
                description="Set the color for the element Fe",
            )
            co = tree.inputs.color(
                "Co",
                (1.0, 0.8509804, 0.5607843, 1.0),
                description="Set the color for the element Co",
            )
            ni = tree.inputs.color(
                "Ni",
                (0.7803922, 0.5411765, 0.5411765, 1.0),
                description="Set the color for the element Ni",
            )
            cu = tree.inputs.color(
                "Cu",
                (0.7843137, 0.5019608, 0.2, 1.0),
                description="Set the color for the element Cu",
            )
            zn = tree.inputs.color(
                "Zn",
                (0.4901961, 0.5019608, 0.6901961, 1.0),
                description="Set the color for the element Zn",
            )
            ga = tree.inputs.color(
                "Ga",
                (0.7607843, 0.5607843, 0.5607843, 1.0),
                description="Set the color for the element Ga",
            )
            ge = tree.inputs.color(
                "Ge",
                (0.4, 0.5607843, 0.5607843, 1.0),
                description="Set the color for the element Ge",
            )
            as_ = tree.inputs.color(
                "As",
                (0.7411765, 0.5019608, 0.890196, 1.0),
                description="Set the color for the element As",
            )
            se = tree.inputs.color(
                "Se",
                (1.0, 0.6313726, 0.0, 1.0),
                description="Set the color for the element Se",
            )
            br = tree.inputs.color(
                "Br",
                (0.6509804, 0.1607843, 0.1607843, 1.0),
                description="Set the color for the element Br",
            )
            kr = tree.inputs.color(
                "Kr",
                (0.3607843, 0.7215686, 0.8196079, 1.0),
                description="Set the color for the element Kr",
            )
            rb = tree.inputs.color(
                "Rb",
                (0.4392157, 0.18039216, 0.6901961, 1.0),
                description="Set the color for the element Rb",
            )
            sr = tree.inputs.color(
                "Sr",
                (0.0, 1.0, 0.0, 1.0),
                description="Set the color for the element Sr",
            )
            y = tree.inputs.color(
                "Y",
                (0.5803922, 1.0, 1.0, 1.0),
                description="Set the color for the element Y",
            )
            zr = tree.inputs.color(
                "Zr",
                (0.5803922, 0.8784314, 0.8784314, 1.0),
                description="Set the color for the element Zr",
            )
        with tree.inputs.panel("41-60", default_closed=True):
            nb = tree.inputs.color(
                "Nb",
                (0.4509804, 0.7607843, 0.7882353, 1.0),
                description="Set the color for the element Nb",
            )
            mo = tree.inputs.color(
                "Mo",
                (0.3294118, 0.709804, 0.709804, 1.0),
                description="Set the color for the element Mo",
            )
            tc = tree.inputs.color(
                "Tc",
                (0.23137255, 0.6196079, 0.6196079, 1.0),
                description="Set the color for the element Tc",
            )
            ru = tree.inputs.color(
                "Ru",
                (0.14117648, 0.4901961, 0.4901961, 1.0),
                description="Set the color for the element Ru",
            )
            rh = tree.inputs.color(
                "Rh",
                (0.03921569, 0.4901961, 0.5490196, 1.0),
                description="Set the color for the element Rh",
            )
            pd = tree.inputs.color(
                "Pd",
                (0.0, 0.4117647, 0.5215687, 1.0),
                description="Set the color for the element Pd",
            )
            ag = tree.inputs.color(
                "Ag",
                (0.7529412, 0.7529412, 0.7529412, 1.0),
                description="Set the color for the element Ag",
            )
            cd = tree.inputs.color(
                "Cd",
                (1.0, 0.8509804, 0.5607843, 1.0),
                description="Set the color for the element Cd",
            )
            in_ = tree.inputs.color(
                "In",
                (0.6509804, 0.4588235, 0.4509804, 1.0),
                description="Set the color for the element In",
            )
            sn = tree.inputs.color(
                "Sn",
                (0.4, 0.5019608, 0.5019608, 1.0),
                description="Set the color for the element Sn",
            )
            sb = tree.inputs.color(
                "Sb",
                (0.6196079, 0.3882353, 0.709804, 1.0),
                description="Set the color for the element Sb",
            )
            te = tree.inputs.color(
                "Te",
                (0.8313726, 0.4784314, 0.0, 1.0),
                description="Set the color for the element Te",
            )
            i = tree.inputs.color(
                "I",
                (0.5803922, 0.0, 0.5803922, 1.0),
                description="Set the color for the element I",
            )
            xe = tree.inputs.color(
                "Xe",
                (0.2588235, 0.6196079, 0.6901961, 1.0),
                description="Set the color for the element Xe",
            )
            cs = tree.inputs.color(
                "Cs",
                (0.3411765, 0.09019608, 0.5607843, 1.0),
                description="Set the color for the element Cs",
            )
            ba = tree.inputs.color(
                "Ba",
                (0.0, 0.7882353, 0.0, 1.0),
                description="Set the color for the element Ba",
            )
            la = tree.inputs.color(
                "La",
                (0.4392157, 0.8313726, 1.0, 1.0),
                description="Set the color for the element La",
            )
            ce = tree.inputs.color(
                "Ce",
                (1.0, 1.0, 0.7803922, 1.0),
                description="Set the color for the element Ce",
            )
            pr = tree.inputs.color(
                "Pr",
                (0.8509804, 1.0, 0.7803922, 1.0),
                description="Set the color for the element Pr",
            )
            nd = tree.inputs.color(
                "Nd",
                (0.7803922, 1.0, 0.7803922, 1.0),
                description="Set the color for the element Nd",
            )
        with tree.inputs.panel("61-80", default_closed=True):
            pm = tree.inputs.color(
                "Pm",
                (0.6392157, 1.0, 0.7803922, 1.0),
                description="Set the color for the element Pm",
            )
            sm = tree.inputs.color(
                "Sm",
                (0.5607843, 1.0, 0.7803922, 1.0),
                description="Set the color for the element Sm",
            )
            eu = tree.inputs.color(
                "Eu",
                (0.3803922, 1.0, 0.7803922, 1.0),
                description="Set the color for the element Eu",
            )
            gd = tree.inputs.color(
                "Gd",
                (0.27058825, 1.0, 0.7803922, 1.0),
                description="Set the color for the element Gd",
            )
            tb = tree.inputs.color(
                "Tb",
                (0.1882353, 1.0, 0.7803922, 1.0),
                description="Set the color for the element Tb",
            )
            dy = tree.inputs.color(
                "Dy",
                (0.12156863, 1.0, 0.7803922, 1.0),
                description="Set the color for the element Dy",
            )
            ho = tree.inputs.color(
                "Ho",
                (0.0, 1.0, 0.6117647, 1.0),
                description="Set the color for the element Ho",
            )
            er = tree.inputs.color(
                "Er",
                (0.0, 0.9019608, 0.4588235, 1.0),
                description="Set the color for the element Er",
            )
            tm = tree.inputs.color(
                "Tm",
                (0.0, 0.8313726, 0.3215686, 1.0),
                description="Set the color for the element Tm",
            )
            yb = tree.inputs.color(
                "Yb",
                (0.0, 0.7490196, 0.21960784, 1.0),
                description="Set the color for the element Yb",
            )
            lu = tree.inputs.color(
                "Lu",
                (0.0, 0.6705883, 0.14117648, 1.0),
                description="Set the color for the element Lu",
            )
            hf = tree.inputs.color(
                "Hf",
                (0.3019608, 0.7607843, 1.0, 1.0),
                description="Set the color for the element Hf",
            )
            ta = tree.inputs.color(
                "Ta",
                (0.3019608, 0.6509804, 1.0, 1.0),
                description="Set the color for the element Ta",
            )
            w = tree.inputs.color(
                "W",
                (0.12941177, 0.5803922, 0.8392157, 1.0),
                description="Set the color for the element W",
            )
            re = tree.inputs.color(
                "Re",
                (0.1490196, 0.4901961, 0.6705883, 1.0),
                description="Set the color for the element Re",
            )
            os = tree.inputs.color(
                "Os",
                (0.1490196, 0.4, 0.5882353, 1.0),
                description="Set the color for the element Os",
            )
            ir = tree.inputs.color(
                "Ir",
                (0.09019608, 0.3294118, 0.5294118, 1.0),
                description="Set the color for the element Ir",
            )
            pt = tree.inputs.color(
                "Pt",
                (0.8156863, 0.8156863, 0.8784314, 1.0),
                description="Set the color for the element Pt",
            )
            au = tree.inputs.color(
                "Au",
                (1.0, 0.8196079, 0.1372549, 1.0),
                description="Set the color for the element Au",
            )
            hg = tree.inputs.color(
                "Hg",
                (0.7215686, 0.7215686, 0.8156863, 1.0),
                description="Set the color for the element Hg",
            )
        color = tree.outputs.color(
            "Color",
            (0.0, 0.0, 0.0, 1.0),
            description="The selected colors based on the `atomic_number`",
        )

        (
            g.IndexSwitch.color(
                AtomicNumber(),
                (
                    (0.8000075, 0.18144128, 0.5491766, 1.0),
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
            >> color
        )


ASSET = ColorElement

ASSET_METADATA = {
    "catalog_id": "d3f975df-8408-4972-a669-8187a57e01d0",
}
