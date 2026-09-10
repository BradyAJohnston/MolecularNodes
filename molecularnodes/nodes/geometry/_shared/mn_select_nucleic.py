# Node group '.MN_select_nucleic' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import BooleanSocket, CustomGeometryGroup, SocketAccessor
from .mn_constants_atom_name_nucleic import MN_constants_atom_name_nucleic


class MN_select_nucleic(CustomGeometryGroup):
    """
    .MN_select_nucleic

    Outputs
    -------
    o.is_backbone : BooleanSocket
        True for atoms that are part of the sugar-phosphate backbone for the nucleotides
    o.is_side_chain : BooleanSocket
        True for atoms that are part of the bases for nucleotides.
    o.is_nucleic : BooleanSocket
        True if the atoms are part of a nucleic acid
    """

    _name = ".MN_select_nucleic"
    _tree_properties = {"node_tool_idname": "geometry._mn_select_nucleic"}

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        is_backbone: BooleanSocket
        """True for atoms that are part of the sugar-phosphate backbone for the nucleotides"""
        is_side_chain: BooleanSocket
        """True for atoms that are part of the bases for nucleotides."""
        is_nucleic: BooleanSocket
        """True if the atoms are part of a nucleic acid"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self):
        super().__init__()

    def _build_group(self, tree):
        is_backbone = tree.outputs.boolean(
            "Is Backbone",
            description="True for atoms that are part of the sugar-phosphate backbone for the nucleotides",
        )
        is_side_chain = tree.outputs.boolean(
            "Is Side Chain",
            description="True for atoms that are part of the bases for nucleotides.",
        )
        is_nucleic = tree.outputs.boolean(
            "Is Nucleic", description="True if the atoms are part of a nucleic acid"
        )

        group = MN_constants_atom_name_nucleic()
        named_attribute = g.NamedAttribute.integer("atom_name")
        (
            (
                (named_attribute.o.attribute >= group.o.backbone_lower)
                & (named_attribute.o.attribute <= group.o.backbone_upper)
            )
            >> is_backbone
        )
        (
            (
                (named_attribute.o.attribute >= group.o.side_chain_lower)
                & (named_attribute.o.attribute <= group.o.side_chain_upper)
            )
            >> is_side_chain
        )
        (
            (
                (named_attribute.o.attribute >= group.o.backbone_lower)
                & (named_attribute.o.attribute <= group.o.side_chain_upper)
            )
            >> is_nucleic
        )
