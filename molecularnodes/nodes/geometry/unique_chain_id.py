# Node-group asset "Unique Chain ID" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputFloat
from .angstrom_to_world import AngstromToWorld
from .chain_id import ChainID
from .integer_run import IntegerRun
from .offset_vector import OffsetVector


class UniqueChainID(AssetGeometryGroup):
    """
    Compute a unique Group ID based on the `chain_id` attribute, but also incrememnting if subsequent points are over a cutoff distance from each other

    Parameters
    ----------
    cutoff : InputFloat
        Threshold distance over which the next points are considered to be part of a new chain

    Inputs
    ------
    i.cutoff : FloatSocket
        Threshold distance over which the next points are considered to be part of a new chain

    Outputs
    -------
    o.group_id : IntegerSocket
        Calculated `Group ID` attribute
    """

    _name = "Unique Chain ID"
    _asset_name = "Unique Chain ID"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {
        "description": "Compute a unique Group ID based on the `chain_id` attribute, but also incrememnting if subsequent points are over a cutoff distance from each other"
    }

    class _Inputs(SocketAccessor):
        cutoff: FloatSocket
        """Threshold distance over which the next points are considered to be part of a new chain"""

    class _Outputs(SocketAccessor):
        group_id: IntegerSocket
        """Calculated `Group ID` attribute"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        cutoff: InputFloat = 4.5,
    ):
        super().__init__(**{"Cutoff": cutoff})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        cutoff = tree.inputs.float(
            "Cutoff",
            4.5,
            description="Threshold distance over which the next points are considered to be part of a new chain",
            min_value=0.0,
            max_value=10_000.0,
            subtype="DISTANCE",
        )
        group_id = tree.outputs.integer(
            "Group ID", description="Calculated `Group ID` attribute"
        )

        position = g.Position()
        compare = OffsetVector(vector=position, offset=-1).o.value.distance(
            position
        ) > AngstromToWorld(angstrom=cutoff)
        (
            IntegerRun(value=ChainID()).o.group_id
            + g.AccumulateField.point.integer(compare).o.leading
            >> group_id
        )


ASSET = UniqueChainID

ASSET_METADATA = {
    "description": "Compute a unique Group ID based on the `chain_id` attribute, but also incrememnting if subsequent points are over a cutoff distance from each other",
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
