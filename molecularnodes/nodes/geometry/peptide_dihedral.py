# Node-group asset 'Peptide Dihedral' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
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
from ._shared.override_index import OverrideIndex
from .accumulate_axis_rotation import AccumulateAxisRotation
from .atom_name import AtomName
from .chain_id import ChainID
from .chain_parameter import ChainParameter
from .is_peptide import IsPeptide
from .is_side_chain import IsSideChain
from .menu_residue_mask import MenuResidueMask
from .residue_mask import ResidueMask


class PeptideDihedral(AssetGeometryGroup):
    """
    Peptide Dihedral

    Parameters
    ----------
    position : InputVector
        Position
    selection : InputBoolean
        The resulting selection must overlap with this input selection
    phi : InputFloat
        Phi
    psi : InputFloat
        Psi

    Inputs
    ------
    i.position : VectorSocket
        Position
    i.selection : BooleanSocket
        The resulting selection must overlap with this input selection
    i.phi : FloatSocket
        Phi
    i.psi : FloatSocket
        Psi

    Outputs
    -------
    o.position : VectorSocket
        Transformed vector
    """

    _name = "Peptide Dihedral"
    _asset_name = "Peptide Dihedral"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"

    class _Inputs(SocketAccessor):
        position: VectorSocket
        """Position"""
        selection: BooleanSocket
        """The resulting selection must overlap with this input selection"""
        phi: FloatSocket
        """Phi"""
        psi: FloatSocket
        """Psi"""

    class _Outputs(SocketAccessor):
        position: VectorSocket
        """Transformed vector"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        position: InputVector = None,
        selection: InputBoolean = True,
        phi: InputFloat = 0.0,
        psi: InputFloat = 0.0,
    ):
        super().__init__(
            **{"Position": position, "Selection": selection, "Phi": phi, "Psi": psi}
        )

    def _build_group(self, tree):
        position = tree.inputs.vector(
            "Position", (0.0, 0.0, 0.0), hide_value=True, default_input="POSITION"
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="The resulting selection must overlap with this input selection",
            hide_value=True,
        )
        phi = tree.inputs.float("Phi", 0.0, subtype="ANGLE")
        psi = tree.inputs.float("Psi", 0.0, subtype="ANGLE")
        position_1 = tree.outputs.vector(
            "Position", description="Transformed vector", subtype="XYZ"
        )

        group = AtomName()
        boolean = g.Boolean(boolean=True)
        group_1 = OverrideIndex(
            selection=IsSideChain(include_ca=False).o.selection,
            index=HydrogenBondingPartner(),
            override=ResidueMask(atom_name=2).o.index,
        )
        group_2 = IsPeptide(and_=selection)
        index_switch = g.IndexSwitch.float(
            group,
            (
                0.0,
                0.0,
                g.Switch.float(ChainParameter().o.residue_index, true=phi),
                psi.point.at(MenuResidueMask(atom_name="CA").o.index),
                0.0,
            ),
        )
        group_3 = AccumulateAxisRotation(
            position=position,
            selection=group_2.o.selection,
            pivot=g.IndexSwitch.boolean(group, (False, boolean, boolean, boolean)),
            angle=index_switch,
            group_id=ChainID(),
            transform_index=group_1,
        )
        group_2.o.selection.switch.vector(position, group_3.o.position) >> position_1


ASSET = PeptideDihedral

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
