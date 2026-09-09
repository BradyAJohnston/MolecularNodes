# Node group '.MN_pivot_peptide' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import BooleanSocket, CustomGeometryGroup, SocketAccessor
from ..atom_name import AtomName
from .mn_select_res_name_peptide import MN_select_res_name_peptide


class MN_pivot_peptide(CustomGeometryGroup):
    """
    .MN_pivot_peptide

    Outputs
    -------
    o.pivot_side_chain : BooleanSocket
        Pivot Side Chain
    """

    _name = ".MN_pivot_peptide"
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        pivot_side_chain: BooleanSocket
        """Pivot Side Chain"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self):
        super().__init__()

    def _build_group(self, tree):
        pivot_side_chain = tree.outputs.boolean("Pivot Side Chain")

        group = MN_select_res_name_peptide(
            ala=True,
            arg=True,
            asn=True,
            asp=True,
            glu=True,
            gln=True,
            his=True,
            leu=True,
            lys=True,
            met=True,
            phe=True,
            trp=True,
            tyr=True,
        )
        group_1 = MN_select_res_name_peptide(
            ala=True,
            arg=True,
            asn=True,
            asp=True,
            cys=True,
            glu=True,
            gln=True,
            his=True,
            ile=True,
            leu=True,
            lys=True,
            met=True,
            phe=True,
            ser=True,
            thr=True,
            trp=True,
            tyr=True,
            val=True,
        )
        group_2 = MN_select_res_name_peptide(
            ala=True,
            arg=True,
            asn=True,
            asp=True,
            cys=True,
            glu=True,
            gln=True,
            his=True,
            ile=True,
            leu=True,
            lys=True,
            met=True,
            phe=True,
            ser=True,
            thr=True,
            trp=True,
            tyr=True,
            val=True,
        )
        group_3 = MN_select_res_name_peptide(arg=True)
        (
            g.IndexSwitch.boolean(
                AtomName(),
                (
                    False,
                    False,
                    group_2.o.selection,
                    False,
                    False,
                    group_1.o.selection,
                    group.o.selection,
                    MN_select_res_name_peptide(ile=True).o.selection,
                    False,
                    False,
                    False,
                    False,
                    MN_select_res_name_peptide(
                        arg=True, gln=True, lys=True
                    ).o.selection,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    MN_select_res_name_peptide(met=True).o.selection,
                    MN_select_res_name_peptide(ile=True, lys=True).o.selection,
                    False,
                    False,
                    False,
                    False,
                    group_3.o.selection,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    False,
                    group_3.o.selection,
                    False,
                    False,
                ),
            )
            >> pivot_side_chain
        )
