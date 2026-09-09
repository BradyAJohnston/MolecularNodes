# Node group '.MN_constants_atom_name_peptide' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import CustomGeometryGroup, IntegerSocket, SocketAccessor


class MN_constants_atom_name_peptide(CustomGeometryGroup):
    """
    .MN_constants_atom_name_peptide

    Outputs
    -------
    o.backbone_lower : IntegerSocket
        Backbone Lower
    o.backbone_upper : IntegerSocket
        Backbone Upper
    o.side_chain_lower : IntegerSocket
        Side Chain Lower
    o.side_chain_upper : IntegerSocket
        Side Chain Upper
    o.alpha_carbon : IntegerSocket
        Alpha Carbon
    """

    _name = ".MN_constants_atom_name_peptide"
    _tree_properties = {"node_tool_idname": "geometry._mn_constants_atom_name_peptide"}

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        backbone_lower: IntegerSocket
        """Backbone Lower"""
        backbone_upper: IntegerSocket
        """Backbone Upper"""
        side_chain_lower: IntegerSocket
        """Side Chain Lower"""
        side_chain_upper: IntegerSocket
        """Side Chain Upper"""
        alpha_carbon: IntegerSocket
        """Alpha Carbon"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self):
        super().__init__()

    def _build_group(self, tree):
        backbone_lower = tree.outputs.integer("Backbone Lower")
        backbone_upper = tree.outputs.integer("Backbone Upper")
        side_chain_lower = tree.outputs.integer("Side Chain Lower")
        side_chain_upper = tree.outputs.integer("Side Chain Upper")
        alpha_carbon = tree.outputs.integer("Alpha Carbon")

        integer = g.Integer(integer=5)
        integer_1 = g.Integer(integer=49)
        integer_2 = g.Integer(integer=4)
        integer_3 = g.Integer()
        integer_4 = g.Integer(integer=2)

        integer_3 >> backbone_lower
        integer_2 >> backbone_upper
        integer >> side_chain_lower
        integer_1 >> side_chain_upper
        integer_4 >> alpha_carbon
