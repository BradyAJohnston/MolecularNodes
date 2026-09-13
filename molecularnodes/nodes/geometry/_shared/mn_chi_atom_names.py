# Node group ".MN_chi_atom_names" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
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

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        x1 = tree.outputs.integer("X1")
        x2 = tree.outputs.integer("X2")
        x3 = tree.outputs.integer("X3")
        x4 = tree.outputs.integer("X4")
        x5 = tree.outputs.integer("X5")

        menu_atom_name = MenuAtomName(atom_name="CG")
        menu_atom_name_1 = MenuAtomName(atom_name="CG1")
        menu_atom_name_2 = MenuAtomName(atom_name="CD")
        menu_atom_name_3 = MenuAtomName(atom_name="OD1")
        menu_atom_name_4 = MenuAtomName(atom_name="CD1")
        menu_atom_name_5 = MenuAtomName(atom_name="OE1")
        menu_atom_name_6 = MenuAtomName(atom_name="CE")
        _menu_residue_name = MenuResidueName()
        _residue_name = ResidueName()
        (
            SwitchResidueName(
                ala=menu_atom_name.o[0],
                arg=menu_atom_name.o[0],
                asn=menu_atom_name.o[0],
                asp=menu_atom_name.o[0],
                cys=MenuAtomName(atom_name="SG").o[0],
                glu=menu_atom_name.o[0],
                gln=menu_atom_name.o[0],
                gly=menu_atom_name.o[0],
                his=menu_atom_name.o[0],
                ile=menu_atom_name_1.o[0],
                leu=menu_atom_name.o[0],
                lys=menu_atom_name.o[0],
                met=menu_atom_name.o[0],
                phe=menu_atom_name.o[0],
                pro=menu_atom_name.o[0],
                ser=MenuAtomName(atom_name="OG").o[0],
                thr=MenuAtomName(atom_name="CG2").o[0],
                trp=menu_atom_name.o[0],
                tyr=menu_atom_name.o[0],
                val=menu_atom_name_1.o[0],
            )
            >> x1
        )
        (
            SwitchResidueName(
                ala=menu_atom_name_2.o[0],
                arg=menu_atom_name_2.o[0],
                asn=menu_atom_name_3.o[0],
                asp=menu_atom_name_3.o[0],
                cys=menu_atom_name_2.o[0],
                glu=menu_atom_name_2.o[0],
                gln=menu_atom_name_2.o[0],
                gly=menu_atom_name_2.o[0],
                his=MenuAtomName(atom_name="CD2").o[0],
                ile=menu_atom_name_2.o[0],
                leu=menu_atom_name_4.o[0],
                lys=menu_atom_name_2.o[0],
                met=MenuAtomName(atom_name="SD").o[0],
                phe=menu_atom_name_4.o[0],
                pro=menu_atom_name_2.o[0],
                ser=menu_atom_name_2.o[0],
                thr=menu_atom_name_2.o[0],
                trp=menu_atom_name_4.o[0],
                tyr=menu_atom_name_4.o[0],
                val=menu_atom_name_2.o[0],
            )
            >> x2
        )
        (
            SwitchResidueName(
                arg=MenuAtomName(atom_name="CZ").o[0],
                lys=MenuAtomName(atom_name="NZ").o[0],
            )
            >> x4
        )
        SwitchResidueName(lys=MenuAtomName(atom_name="NH1").o[0]) >> x5
        (
            SwitchResidueName(
                ala=0,
                arg=MenuAtomName(atom_name="NE").o[0],
                asn=0,
                asp=0,
                cys=0,
                glu=menu_atom_name_5.o[0],
                gln=menu_atom_name_5.o[0],
                gly=0,
                his=0,
                ile=0,
                leu=0,
                lys=menu_atom_name_6.o[0],
                met=menu_atom_name_6.o[0],
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
