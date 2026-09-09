# Node group '.MN_chi_atom_names' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy.builder import CustomGeometryGroup, IntegerSocket, SocketAccessor
from ..menu_atom_name import MenuAtomName
from ..menu_residue_name import MenuResidueName
from ..residue_name import ResidueName
from ..switch_residue_name import SwitchResidueName


class MN_chi_atom_names(CustomGeometryGroup):
    """
    .MN_chi_atom_names

    Outputs
    -------
    o.x1 : IntegerSocket
        X1
    o.x2 : IntegerSocket
        X2
    o.x3 : IntegerSocket
        X3
    o.x4 : IntegerSocket
        X4
    o.x5 : IntegerSocket
        X5
    """

    _name = ".MN_chi_atom_names"

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        x1: IntegerSocket
        """X1"""
        x2: IntegerSocket
        """X2"""
        x3: IntegerSocket
        """X3"""
        x4: IntegerSocket
        """X4"""
        x5: IntegerSocket
        """X5"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self):
        super().__init__()

    def _build_group(self, tree):
        x1 = tree.outputs.integer("X1")
        x2 = tree.outputs.integer("X2")
        x3 = tree.outputs.integer("X3")
        x4 = tree.outputs.integer("X4")
        x5 = tree.outputs.integer("X5")

        group = MenuAtomName(atom_name="CG")
        group_1 = MenuAtomName(atom_name="CG1")
        (
            SwitchResidueName(
                ala=group.o.socket_1,
                arg=group.o.socket_1,
                asn=group.o.socket_1,
                asp=group.o.socket_1,
                cys=MenuAtomName(atom_name="SG").o.socket_1,
                glu=group.o.socket_1,
                gln=group.o.socket_1,
                gly=group.o.socket_1,
                his=group.o.socket_1,
                ile=group_1.o.socket_1,
                leu=group.o.socket_1,
                lys=group.o.socket_1,
                met=group.o.socket_1,
                phe=group.o.socket_1,
                pro=group.o.socket_1,
                ser=MenuAtomName(atom_name="OG").o.socket_1,
                thr=MenuAtomName(atom_name="CG2").o.socket_1,
                trp=group.o.socket_1,
                tyr=group.o.socket_1,
                val=group_1.o.socket_1,
            )
            >> x1
        )
        group_2 = MenuAtomName(atom_name="CD")
        group_3 = MenuAtomName(atom_name="OD1")
        group_4 = MenuAtomName(atom_name="CD1")
        (
            SwitchResidueName(
                ala=group_2.o.socket_1,
                arg=group_2.o.socket_1,
                asn=group_3.o.socket_1,
                asp=group_3.o.socket_1,
                cys=group_2.o.socket_1,
                glu=group_2.o.socket_1,
                gln=group_2.o.socket_1,
                gly=group_2.o.socket_1,
                his=MenuAtomName(atom_name="CD2").o.socket_1,
                ile=group_2.o.socket_1,
                leu=group_4.o.socket_1,
                lys=group_2.o.socket_1,
                met=MenuAtomName(atom_name="SD").o.socket_1,
                phe=group_4.o.socket_1,
                pro=group_2.o.socket_1,
                ser=group_2.o.socket_1,
                thr=group_2.o.socket_1,
                trp=group_4.o.socket_1,
                tyr=group_4.o.socket_1,
                val=group_2.o.socket_1,
            )
            >> x2
        )
        _group_5 = ResidueName()
        (
            SwitchResidueName(
                arg=MenuAtomName(atom_name="CZ").o.socket_1,
                lys=MenuAtomName(atom_name="NZ").o.socket_1,
            )
            >> x4
        )
        SwitchResidueName(lys=MenuAtomName(atom_name="NH1").o.socket_1) >> x5
        group_6 = MenuAtomName(atom_name="OE1")
        group_7 = MenuAtomName(atom_name="CE")
        (
            SwitchResidueName(
                ala=0,
                arg=MenuAtomName(atom_name="NE").o.socket_1,
                asn=0,
                asp=0,
                cys=0,
                glu=group_6.o.socket_1,
                gln=group_6.o.socket_1,
                gly=0,
                his=0,
                ile=0,
                leu=0,
                lys=group_7.o.socket_1,
                met=group_7.o.socket_1,
                phe=0,
                pro=0,
                ser=0,
                thr=0,
                trp=0,
                tyr=0,
                val=0,
            )
            >> x3
        )
        _group_8 = MenuResidueName()
