# Node-group asset "Peptide Chi" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputBoolean, InputFloat, InputVector
from ._shared.hydrogen_bonding_partner import HydrogenBondingPartner
from ._shared.mn_peptide_chi_values import MN_peptide_chi_values
from ._shared.mn_pivot_peptide import MN_pivot_peptide
from ._shared.override_index import OverrideIndex
from .accumulate_axis_rotation import AccumulateAxisRotation
from .is_peptide import IsPeptide
from .menu_atom_name import MenuAtomName
from .menu_residue_mask import MenuResidueMask
from .ures_id import UResID


class PeptideChi(AssetGeometryGroup):
    """
    Peptide Chi

    Parameters
    ----------
    position : InputVector
        Position
    selection : InputBoolean
        Selection
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
    i.position : VectorSocket
        Position
    i.selection : BooleanSocket
        Selection
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
    o.position : VectorSocket
        Position
    """

    _name = "Peptide Chi"
    _asset_name = "Peptide Chi"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"

    class _Inputs(SocketAccessor):
        position: VectorSocket
        """Position"""
        selection: BooleanSocket
        """Selection"""
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
        position: VectorSocket
        """Position"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        position: InputVector = None,
        selection: InputBoolean = True,
        x1: InputFloat = 0.0,
        x2: InputFloat = 0.0,
        x3: InputFloat = 0.0,
        x4: InputFloat = 0.0,
        x5: InputFloat = 0.0,
    ):
        super().__init__(
            **{
                "Position": position,
                "Selection": selection,
                "X1": x1,
                "X2": x2,
                "X3": x3,
                "X4": x4,
                "X5": x5,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        position = tree.inputs.vector(
            "Position", (0.0, 0.0, 0.0), subtype="XYZ", default_input="POSITION"
        )
        selection = tree.inputs.boolean("Selection", True, hide_value=True)
        x1 = tree.inputs.float("X1", 0.0, subtype="ANGLE")
        x2 = tree.inputs.float("X2", 0.0, subtype="ANGLE")
        x3 = tree.inputs.float("X3", 0.0, subtype="ANGLE")
        x4 = tree.inputs.float("X4", 0.0, subtype="ANGLE")
        x5 = tree.inputs.float("X5", 0.0, subtype="ANGLE")
        position_1 = tree.outputs.vector("Position", subtype="XYZ")

        group = MN_pivot_peptide()
        group_1 = OverrideIndex(
            selection=MenuAtomName(atom_name="CG2").o.selection,
            override=MenuResidueMask(atom_name="CB").o.index,
        )
        group_2 = AccumulateAxisRotation(
            position=position,
            selection=IsPeptide(and_=selection).o.selection & group,
            pivot=group,
            angle=MN_peptide_chi_values(x1=x1, x2=x2, x3=x3, x4=x4, x5=x5).o.value,
            group_id=UResID().o.ures_id,
            transform_index=group_1.o.output.point.at(HydrogenBondingPartner()),
        )
        (
            (~MenuAtomName(atom_name="OXT").o.selection & selection).switch.vector(
                position, group_2.o.position
            )
            >> position_1
        )


ASSET = PeptideChi

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
