# Node group ".MN_peptide_chi_values" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    CustomGeometryGroup,
    FloatSocket,
    IntegerSocket,
    SocketAccessor,
)
from nodebpy.types import InputFloat
from ..atom_name import AtomName
from ..menu_residue_name import MenuResidueName


class MN_peptide_chi_values(CustomGeometryGroup):
    """
    .MN_peptide_chi_values

    Parameters
    ----------
    x1 : InputFloat
        X1
    x2 : InputFloat
        X2
    x3 : InputFloat
        X3
    x4 : InputFloat
        X4
    x5 : InputFloat
        X5

    Inputs
    ------
    i.x1 : FloatSocket
        X1
    i.x2 : FloatSocket
        X2
    i.x3 : FloatSocket
        X3
    i.x4 : FloatSocket
        X4
    i.x5 : FloatSocket
        X5

    Outputs
    -------
    o.value : FloatSocket
        Value
    o.atom_name : IntegerSocket
        atom_name
    """

    _name = ".MN_peptide_chi_values"

    class _Inputs(SocketAccessor):
        x1: FloatSocket
        """X1"""
        x2: FloatSocket
        """X2"""
        x3: FloatSocket
        """X3"""
        x4: FloatSocket
        """X4"""
        x5: FloatSocket
        """X5"""

    class _Outputs(SocketAccessor):
        value: FloatSocket
        """Value"""
        atom_name: IntegerSocket

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        x1: InputFloat = 0.0,
        x2: InputFloat = 0.0,
        x3: InputFloat = 0.0,
        x4: InputFloat = 0.0,
        x5: InputFloat = 0.0,
    ):
        super().__init__(**{"X1": x1, "X2": x2, "X3": x3, "X4": x4, "X5": x5})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        x1 = tree.inputs.float("X1", 0.0)
        x2 = tree.inputs.float("X2", 0.0)
        x3 = tree.inputs.float("X3", 0.0)
        x4 = tree.inputs.float("X4", 0.0)
        x5 = tree.inputs.float("X5", 0.0)
        value = tree.outputs.float("Value")
        atom_name = tree.outputs.integer("atom_name")

        group = AtomName()
        (
            g.IndexSwitch.integer(
                group,
                (
                    0,
                    0,
                    0,
                    0,
                    0,
                    2,
                    5,
                    5,
                    0,
                    0,
                    0,
                    0,
                    6,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    12,
                    0,
                    0,
                    0,
                    0,
                    12,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    25,
                ),
            )
            >> atom_name
        )
        index_switch = g.IndexSwitch.float(
            group,
            (
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                x1,
                x2,
                x2,
                0.0,
                0.0,
                0.0,
                0.0,
                x3,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                MenuResidueName(residue_name="MET").o.selection.switch.float(x4, x3),
                x4,
                0.0,
                0.0,
                0.0,
                0.0,
                x4,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                MenuResidueName(residue_name="TYR").o.selection.switch.float(x5),
            ),
        )
        (
            MenuResidueName(residue_name="PRO").o.selection.switch.float(index_switch)
            >> value
        )
