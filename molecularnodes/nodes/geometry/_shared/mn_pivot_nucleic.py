# Node group '.MN_pivot_nucleic' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import BooleanSocket, CustomGeometryGroup, SocketAccessor
from ..atom_name import AtomName
from ..select_nucleic_type import SelectNucleicType


class MN_pivot_nucleic(CustomGeometryGroup):
    """
    .MN_pivot_nucleic

    Outputs
    -------
    o.accumulate_backbone : BooleanSocket
        Accumulate Backbone
    o.pivot_backbone : BooleanSocket
        Pivot Backbone
    o.accumulate_base : BooleanSocket
        Accumulate Base
    o.pivot_base : BooleanSocket
        Pivot Base
    """

    _name = ".MN_pivot_nucleic"

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        accumulate_backbone: BooleanSocket
        """Accumulate Backbone"""
        pivot_backbone: BooleanSocket
        """Pivot Backbone"""
        accumulate_base: BooleanSocket
        """Accumulate Base"""
        pivot_base: BooleanSocket
        """Pivot Base"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self):
        super().__init__()

    def _build_group(self, tree):
        accumulate_backbone = tree.outputs.boolean("Accumulate Backbone")
        pivot_backbone = tree.outputs.boolean("Pivot Backbone")
        accumulate_base = tree.outputs.boolean("Accumulate Base")
        pivot_base = tree.outputs.boolean("Pivot Base")

        group = AtomName()
        boolean_math = (
            g.Compare.integer.equal(group, 50).o.result
            | g.Compare.integer.equal(group, 53)
            | g.Compare.integer.equal(group, 54)
        )
        boolean_math_1 = (
            boolean_math
            | g.Compare.integer.equal(group, 55)
            | g.Compare.integer.equal(group, 58)
        )
        (boolean_math_1 | g.Compare.integer.equal(group, 57)) >> pivot_backbone
        group_1 = AtomName()
        g.Compare.integer.equal(group_1, 61) >> pivot_base
        (
            g.Compare.integer.equal(
                SelectNucleicType().o.is_pyrimidine.switch.integer(63, 62), group_1
            )
            >> accumulate_base
        )

        boolean_math_1 >> accumulate_backbone
