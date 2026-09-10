# Node-group asset "Nucleic Dihedral" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
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
from ._shared.mn_pivot_nucleic import MN_pivot_nucleic
from ._shared.override_index import OverrideIndex
from .accumulate_axis_rotation import AccumulateAxisRotation
from .atom_name import AtomName
from .bond_count import BondCount
from .chain_id import ChainID
from .chain_parameter import ChainParameter
from .is_nucleic import IsNucleic
from .residue_mask import ResidueMask


class NucleicDihedral(AssetGeometryGroup):
    """
    Nucleic Dihedral

    Parameters
    ----------
    position : InputVector
        Position
    selection : InputBoolean
        The resulting selection must overlap with this input selection
    alpha : InputFloat
        Alpha
    beta : InputFloat
        Beta
    gamma : InputFloat
        Gamma
    epsilon : InputFloat
        Epsilon
    zeta : InputFloat
        Amount to rotate around the axis

    Inputs
    ------
    i.position : VectorSocket
        Position
    i.selection : BooleanSocket
        The resulting selection must overlap with this input selection
    i.alpha : FloatSocket
        Alpha
    i.beta : FloatSocket
        Beta
    i.gamma : FloatSocket
        Gamma
    i.epsilon : FloatSocket
        Epsilon
    i.zeta : FloatSocket
        Amount to rotate around the axis

    Outputs
    -------
    o.position : VectorSocket
        Transformed vector
    """

    _name = "Nucleic Dihedral"
    _asset_name = "Nucleic Dihedral"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"

    class _Inputs(SocketAccessor):
        position: VectorSocket
        """Position"""
        selection: BooleanSocket
        """The resulting selection must overlap with this input selection"""
        alpha: FloatSocket
        """Alpha"""
        beta: FloatSocket
        """Beta"""
        gamma: FloatSocket
        """Gamma"""
        epsilon: FloatSocket
        """Epsilon"""
        zeta: FloatSocket
        """Amount to rotate around the axis"""

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
        alpha: InputFloat = 0.0,
        beta: InputFloat = 0.0,
        gamma: InputFloat = 0.0,
        epsilon: InputFloat = 0.0,
        zeta: InputFloat = 0.0,
    ):
        super().__init__(
            **{
                "Position": position,
                "Selection": selection,
                "Alpha": alpha,
                "Beta": beta,
                "Gamma": gamma,
                "Epsilon": epsilon,
                "Zeta": zeta,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        position = tree.inputs.vector(
            "Position", (0.0, 0.0, 0.0), hide_value=True, default_input="POSITION"
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="The resulting selection must overlap with this input selection",
            hide_value=True,
        )
        alpha = tree.inputs.float("Alpha", 0.0, subtype="ANGLE")
        beta = tree.inputs.float("Beta", 0.0, subtype="ANGLE")
        gamma = tree.inputs.float("Gamma", 0.0, subtype="ANGLE")
        epsilon = tree.inputs.float("Epsilon", 0.0, subtype="ANGLE")
        zeta = tree.inputs.float(
            "Zeta", 0.0, description="Amount to rotate around the axis", subtype="ANGLE"
        )
        position_1 = tree.outputs.vector(
            "Position", description="Transformed vector", subtype="XYZ"
        )

        group = IsNucleic(and_=selection)
        group_1 = AtomName()
        group_2 = OverrideIndex(
            selection=(group_1 > 58) & (group_1 <= 115)
            | g.Compare.integer.equal(group_1, 56),
            override=ResidueMask(atom_name=55).o.index,
        )
        index_switch = g.IndexSwitch.float(
            group_1.o.atom_name - 50,
            (
                zeta,
                0.0,
                0.0,
                alpha,
                g.Switch.float(ChainParameter().o.residue_index, true=beta),
                gamma,
                0.0,
                0.0,
                epsilon,
                0.0,
                0.0,
            ),
        )
        group_3 = AccumulateAxisRotation(
            position=position,
            selection=group.o.selection,
            pivot=MN_pivot_nucleic().o.pivot_backbone,
            angle=(BondCount().o.bonds > 1).switch.float(true=index_switch),
            group_id=ChainID(),
            transform_index=group_2,
        )
        group.o.selection.switch.vector(position, group_3.o.position) >> position_1


ASSET = NucleicDihedral

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
