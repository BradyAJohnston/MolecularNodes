# Node-group asset "Chain Parameter" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from ._shared.index_to_factor import IndexToFactor
from .chain_id import ChainID
from .group_info import GroupInfo
from .sub_group_info import SubGroupInfo
from .ures_id import UResID


class ChainParameter(AssetGeometryGroup):
    """
    Information for each residue within the context of the chain

    Outputs
    -------
    o.factor : FloatSocket
        A residue's relative position along its chain, 0 for the first residue and 1 for the last
    o.residue_index : IntegerSocket
        Index of the residue within its chain, starting at 0 for the first residue
    o.residue_count : IntegerSocket
        Number of residues in the chain
    o.is_first : BooleanSocket
        If the atom is in the first residue of its chain
    o.is_last : BooleanSocket
        If the atom is in the last residue of its chain
    o.index_of_first : IntegerSocket
        `Index` (in the whole structure) of the first atom in the chain
    o.index_of_last : IntegerSocket
        `Index` (in the whole structure) of the last atom in the chain
    """

    _name = "Chain Parameter"
    _asset_name = "Chain Parameter"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"
    _tree_properties = {
        "description": "Information for each residue within the context of the chain",
        "node_tool_idname": "geometry.chain_parameter",
    }

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        factor: FloatSocket
        """A residue's relative position along its chain, 0 for the first residue and 1 for the last"""
        residue_index: IntegerSocket
        """Index of the residue within its chain, starting at 0 for the first residue"""
        residue_count: IntegerSocket
        """Number of residues in the chain"""
        is_first: BooleanSocket
        """If the atom is in the first residue of its chain"""
        is_last: BooleanSocket
        """If the atom is in the last residue of its chain"""
        index_of_first: IntegerSocket
        """`Index` (in the whole structure) of the first atom in the chain"""
        index_of_last: IntegerSocket
        """`Index` (in the whole structure) of the last atom in the chain"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self) -> None:
        super().__init__()

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        factor = tree.outputs.float(
            "Factor",
            description="A residue's relative position along its chain, 0 for the first residue and 1 for the last",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        residue_index = tree.outputs.integer(
            "Residue Index",
            description="Index of the residue within its chain, starting at 0 for the first residue",
        )
        residue_count = tree.outputs.integer(
            "Residue Count", description="Number of residues in the chain"
        )
        is_first = tree.outputs.boolean(
            "Is First", description="If the atom is in the first residue of its chain"
        )
        is_last = tree.outputs.boolean(
            "Is Last", description="If the atom is in the last residue of its chain"
        )
        index_of_first = tree.outputs.integer(
            "Index of First",
            description="`Index` (in the whole structure) of the first atom in the chain",
        )
        index_of_last = tree.outputs.integer(
            "Index of Last",
            description="`Index` (in the whole structure) of the last atom in the chain",
        )

        chain_id = ChainID()
        group_info = GroupInfo(group_id=chain_id)
        sub_group_info = SubGroupInfo(
            sub_group_id=UResID().o.ures_id, group_id=chain_id
        )
        g.Compare.integer.equal(sub_group_info.o.sub_group_id, 0) >> is_first
        (
            IndexToFactor(
                index=sub_group_info.o.sub_group_id,
                size=sub_group_info.o.sub_group_total,
            )
            >> factor
        )
        (
            g.Compare.integer.equal(
                sub_group_info.o.sub_group_id, sub_group_info.o.sub_group_total - 1
            )
            >> is_last
        )

        sub_group_info.o.sub_group_id >> residue_index
        sub_group_info.o.sub_group_total >> residue_count
        group_info.o.index_of_first >> index_of_first
        group_info.o.index_of_last >> index_of_last


ASSET = ChainParameter

ASSET_METADATA = {
    "description": "Information for each residue within the context of the chain",
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
