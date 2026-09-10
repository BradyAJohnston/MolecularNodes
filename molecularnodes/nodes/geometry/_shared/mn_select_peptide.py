# Node group '.MN_select_peptide' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import BooleanSocket, CustomGeometryGroup, SocketAccessor
from ..atom_name import AtomName
from .mn_constants_atom_name_peptide import MN_constants_atom_name_peptide


class MN_select_peptide(CustomGeometryGroup):
    """
    .MN_select_peptide

    Outputs
    -------
    o.is_backbone : BooleanSocket
        Is Backbone
    o.is_side_chain : BooleanSocket
        Is Side Chain
    o.is_peptide : BooleanSocket
        Is Peptide
    o.is_alpha_carbon : BooleanSocket
        Is Alpha Carbon
    """

    _name = ".MN_select_peptide"
    _tree_properties = {"node_tool_idname": "geometry._mn_select_peptide"}

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        is_backbone: BooleanSocket
        """Is Backbone"""
        is_side_chain: BooleanSocket
        """Is Side Chain"""
        is_peptide: BooleanSocket
        """Is Peptide"""
        is_alpha_carbon: BooleanSocket
        """Is Alpha Carbon"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self):
        super().__init__()

    def _build_group(self, tree):
        is_backbone = tree.outputs.boolean("Is Backbone")
        is_side_chain = tree.outputs.boolean("Is Side Chain")
        is_peptide = tree.outputs.boolean("Is Peptide")
        is_alpha_carbon = tree.outputs.boolean("Is Alpha Carbon")

        group = MN_constants_atom_name_peptide()
        group_1 = AtomName()
        (
            ((group_1 >= group.o.backbone_lower) & (group_1 <= group.o.backbone_upper))
            >> is_backbone
        )
        (
            (
                (group_1 >= group.o.backbone_lower)
                & (group_1 <= group.o.side_chain_upper)
            )
            >> is_peptide
        )
        compare = g.Compare.integer.equal(group_1, group.o.alpha_carbon)
        (
            (
                (group_1 >= group.o.side_chain_lower)
                & (group_1 <= group.o.side_chain_upper)
                | compare
            )
            >> is_side_chain
        )

        compare >> is_alpha_carbon
