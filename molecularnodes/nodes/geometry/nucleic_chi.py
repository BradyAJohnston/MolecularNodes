# Node-group asset 'Nucleic Chi' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
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
from .accumulate_axis_rotation import AccumulateAxisRotation
from .is_nucleic import IsNucleic
from .unique_residue_id import UniqueResidueID


class NucleicChi(AssetGeometryGroup):
    """
    Nucleic Chi

    Parameters
    ----------
    position : InputVector
        Position
    selection : InputBoolean
        The resulting selection must overlap with this input selection
    x1 : InputFloat
        Amount to rotate around the axis

    Inputs
    ------
    i.position : VectorSocket
        Position
    i.selection : BooleanSocket
        The resulting selection must overlap with this input selection
    i.x1 : FloatSocket
        Amount to rotate around the axis

    Outputs
    -------
    o.position : VectorSocket
        Transformed vector
    """

    _name = "Nucleic Chi"
    _asset_name = "Nucleic Chi"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"

    class _Inputs(SocketAccessor):
        position: VectorSocket
        """Position"""
        selection: BooleanSocket
        """The resulting selection must overlap with this input selection"""
        x1: FloatSocket
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
        x1: InputFloat = 0.0,
    ):
        super().__init__(**{"Position": position, "Selection": selection, "X1": x1})

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
        x1 = tree.inputs.float(
            "X1", 0.0, description="Amount to rotate around the axis", subtype="ANGLE"
        )
        position_1 = tree.outputs.vector(
            "Position", description="Transformed vector", subtype="XYZ"
        )

        group = IsNucleic(and_=selection)
        group_1 = MN_pivot_nucleic()
        group_2 = AccumulateAxisRotation(
            position=position,
            selection=group.o.selection,
            pivot=group_1.o.pivot_base,
            angle=group_1.o.accumulate_base.switch.float(true=x1),
            group_id=UniqueResidueID(),
        )
        group.o.selection.switch.vector(position, group_2.o.position) >> position_1


ASSET = NucleicChi

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
