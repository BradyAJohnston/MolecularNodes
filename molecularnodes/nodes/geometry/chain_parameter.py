# Node-group asset 'Chain Parameter' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from ._shared.index_to_factor import IndexToFactor
from .chain_id import ChainID
from .group_info import GroupInfo
from .residue_id import ResidueID
from .sub_group_info import SubGroupInfo


class ChainParameter(AssetGeometryGroup):
    """
    Information for each residue within the context of the chain

    Outputs
    -------
    o.factor : FloatSocket
        A residues relative position along a chain. 0 being the first residue in a chain, 1 being the last
    o.residue_count : IntegerSocket
        Number of residues in the chain
    o.residue_index : IntegerSocket
        Res ID along the chain if counting from 1
    o.first_res_id : IntegerSocket
        The first Res ID in a chain (truncated chains start above 1)
    o.last_res_id : IntegerSocket
        The Res ID of the last residue in chain (not equal to Length if chain is truncated)
    o.index_of_first : IntegerSocket
        Index in whole structure of the first atom in the chain
    o.index_of_last : IntegerSocket
        Index in the whole structure the last atom in the chain
    """

    _name = "Chain Parameter"
    _asset_name = "Chain Parameter"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {
        "description": "Information for each residue within the context of the chain",
        "node_tool_idname": "geometry.chain_parameter",
    }

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        factor: FloatSocket
        """A residues relative position along a chain. 0 being the first residue in a chain, 1 being the last"""
        residue_count: IntegerSocket
        """Number of residues in the chain"""
        residue_index: IntegerSocket
        """Res ID along the chain if counting from 1"""
        first_res_id: IntegerSocket
        """The first Res ID in a chain (truncated chains start above 1)"""
        last_res_id: IntegerSocket
        """The Res ID of the last residue in chain (not equal to Length if chain is truncated)"""
        index_of_first: IntegerSocket
        """Index in whole structure of the first atom in the chain"""
        index_of_last: IntegerSocket
        """Index in the whole structure the last atom in the chain"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self):
        super().__init__()

    def _build_group(self, tree):
        factor = tree.outputs.float(
            "Factor",
            description="A residues relative position along a chain. 0 being the first residue in a chain, 1 being the last",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        residue_count = tree.outputs.integer(
            "Residue Count", description="Number of residues in the chain"
        )
        residue_index = tree.outputs.integer(
            "Residue Index", description="Res ID along the chain if counting from 1"
        )
        first_res_id = tree.outputs.integer(
            "First res_id",
            description="The first Res ID in a chain (truncated chains start above 1)",
        )
        last_res_id = tree.outputs.integer(
            "Last res_id",
            description="The Res ID of the last residue in chain (not equal to Length if chain is truncated)",
        )
        index_of_first = tree.outputs.integer(
            "Index of First",
            description="Index in whole structure of the first atom in the chain",
        )
        index_of_last = tree.outputs.integer(
            "Index of Last",
            description="Index in the whole structure the last atom in the chain",
        )

        group = ChainID()
        group_1 = GroupInfo(group_id=group)
        group_2 = SubGroupInfo(sub_group_id=ResidueID(), group_id=group)
        (
            IndexToFactor(index=group_2.o.sub_group_id, size=group_2.o.sub_group_total)
            >> factor
        )
        ResidueID(index=group_2.o.index_of_first) >> first_res_id
        ResidueID(index=group_2.o.index_of_last) >> last_res_id

        group_2.o.sub_group_total >> residue_count
        group_2.o.sub_group_id >> residue_index
        group_1.o.index_of_first >> index_of_first
        group_1.o.index_of_last >> index_of_last


ASSET = ChainParameter

ASSET_METADATA = {
    "description": "Information for each residue within the context of the chain",
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
