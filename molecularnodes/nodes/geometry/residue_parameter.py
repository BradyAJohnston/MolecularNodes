# Node-group asset "Residue Parameter" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from .sub_group_info import SubGroupInfo
from .ures_id import UResID


class ResidueParameter(AssetGeometryGroup):
    """
    Residue Parameter

    Outputs
    -------
    o.factor : FloatSocket
        An atom's relative position in its residue, 0 for the first atom and 1 for the last
    o.atom_index : IntegerSocket
        Index of the atom within its residue, starting at 0 for the first atom
    o.atom_count : IntegerSocket
        Number of atoms in the residue
    o.is_first : BooleanSocket
        If the atom is the first atom in its residue
    o.is_last : BooleanSocket
        If the atom is the last atom in its residue
    o.index_of_first : IntegerSocket
        `Index` (in the whole structure) of the first atom in the residue
    o.index_of_last : IntegerSocket
        `Index` (in the whole structure) of the last atom in the residue
    """

    _name = "Residue Parameter"
    _asset_name = "Residue Parameter"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.residue_parameter"}

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        factor: FloatSocket
        """An atom's relative position in its residue, 0 for the first atom and 1 for the last"""
        atom_index: IntegerSocket
        """Index of the atom within its residue, starting at 0 for the first atom"""
        atom_count: IntegerSocket
        """Number of atoms in the residue"""
        is_first: BooleanSocket
        """If the atom is the first atom in its residue"""
        is_last: BooleanSocket
        """If the atom is the last atom in its residue"""
        index_of_first: IntegerSocket
        """`Index` (in the whole structure) of the first atom in the residue"""
        index_of_last: IntegerSocket
        """`Index` (in the whole structure) of the last atom in the residue"""

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
            description="An atom's relative position in its residue, 0 for the first atom and 1 for the last",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        atom_index = tree.outputs.integer(
            "Atom Index",
            description="Index of the atom within its residue, starting at 0 for the first atom",
        )
        atom_count = tree.outputs.integer(
            "Atom Count", description="Number of atoms in the residue"
        )
        is_first = tree.outputs.boolean(
            "Is First", description="If the atom is the first atom in its residue"
        )
        is_last = tree.outputs.boolean(
            "Is Last", description="If the atom is the last atom in its residue"
        )
        index_of_first = tree.outputs.integer(
            "Index of First",
            description="`Index` (in the whole structure) of the first atom in the residue",
        )
        index_of_last = tree.outputs.integer(
            "Index of Last",
            description="`Index` (in the whole structure) of the last atom in the residue",
        )

        sub_group_info = SubGroupInfo(
            sub_group_id=UResID().o.ures_id, group_id=ChainID()
        )
        g.Compare.integer.equal(sub_group_info.o.index_in_group_id, 0) >> is_first
        (
            IndexToFactor(
                index=sub_group_info.o.index_in_group_id, size=sub_group_info.o.size
            )
            >> factor
        )
        (
            g.Compare.integer.equal(
                sub_group_info.o.index_in_group_id, sub_group_info.o.size - 1
            )
            >> is_last
        )

        sub_group_info.o.index_in_group_id >> atom_index
        sub_group_info >> atom_count
        sub_group_info.o.index_of_first >> index_of_first
        sub_group_info.o.index_of_last >> index_of_last


ASSET = ResidueParameter

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
