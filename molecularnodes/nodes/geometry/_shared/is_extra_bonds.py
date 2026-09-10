# Node group '.Is Extra Bonds' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import BooleanSocket, CustomGeometryGroup, SocketAccessor


class IsExtraBonds(CustomGeometryGroup):
    """
    .Is Extra Bonds

    Outputs
    -------
    o.is_extra_bonds : BooleanSocket
        Is Extra Bonds
    o.is_double_bond : BooleanSocket
        Is Double Bond
    o.is_triple_bond : BooleanSocket
        Is Triple Bond
    """

    _name = ".Is Extra Bonds"
    _tree_properties = {"node_tool_idname": "geometry._is_extra_bonds"}

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        is_extra_bonds: BooleanSocket
        """Is Extra Bonds"""
        is_double_bond: BooleanSocket
        """Is Double Bond"""
        is_triple_bond: BooleanSocket
        """Is Triple Bond"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self):
        super().__init__()

    def _build_group(self, tree):
        is_extra_bonds = tree.outputs.boolean("Is Extra Bonds")
        is_double_bond = tree.outputs.boolean("Is Double Bond")
        is_triple_bond = tree.outputs.boolean("Is Triple Bond")

        named_attribute = g.NamedAttribute.integer("bond_type")
        compare = g.Compare.integer.equal(named_attribute.o.attribute, 2)
        compare_1 = g.Compare.integer.equal(named_attribute.o.attribute, 3)
        (
            (
                compare.o.result
                | compare_1
                | g.Compare.integer.equal(named_attribute.o.attribute, 6)
            )
            >> is_extra_bonds
        )

        compare >> is_double_bond
        compare_1 >> is_triple_bond
