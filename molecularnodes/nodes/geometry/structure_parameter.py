# Node-group asset "Structure Parameter" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from ._shared.index_to_factor import IndexToFactor
from .group_info import GroupInfo
from .sub_group_info import SubGroupInfo
from .unique_residue_id import UniqueResidueID


class StructureParameter(AssetGeometryGroup):
    """
    Structure Parameter

    Outputs
    -------
    o.atom_factor : FloatSocket
        Factor of the atom within the structure, which is the relative position of the `Index` within the overall structure between 0 and 1
    o.index : IntegerSocket
        `Index` of the point within the structure. Equal to the `Index` input node
    o.atom_count : IntegerSocket
        Number of atoms in the structure, equal to the size of the Point Domain
    o.index_of_first : IntegerSocket
        `Index` of first atom in the entire structure, which will always be 0
    o.index_of_last : IntegerSocket
        Index of last atom in the entire structure. Equal to `Atom Count` - 1
    o.residue_factor : FloatSocket
        Residue Factor
    o.residue_index : IntegerSocket
        A unique `Group ID` that increases whenever `Sub Group ID` or `Group ID` change
    o.residue_count : IntegerSocket
        Residue Count
    """

    _name = "Structure Parameter"
    _asset_name = "Structure Parameter"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        atom_factor: FloatSocket
        """Factor of the atom within the structure, which is the relative position of the `Index` within the overall structure between 0 and 1"""
        index: IntegerSocket
        """`Index` of the point within the structure. Equal to the `Index` input node"""
        atom_count: IntegerSocket
        """Number of atoms in the structure, equal to the size of the Point Domain"""
        index_of_first: IntegerSocket
        """`Index` of first atom in the entire structure, which will always be 0"""
        index_of_last: IntegerSocket
        """Index of last atom in the entire structure. Equal to `Atom Count` - 1"""
        residue_factor: FloatSocket
        """Residue Factor"""
        residue_index: IntegerSocket
        """A unique `Group ID` that increases whenever `Sub Group ID` or `Group ID` change"""
        residue_count: IntegerSocket
        """Residue Count"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self):
        super().__init__()

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atom_factor = tree.outputs.float(
            "Atom Factor",
            description="Factor of the atom within the structure, which is the relative position of the `Index` within the overall structure between 0 and 1",
            subtype="FACTOR",
        )
        index = tree.outputs.integer(
            "Index",
            description="`Index` of the point within the structure. Equal to the `Index` input node",
        )
        atom_count = tree.outputs.integer(
            "Atom Count",
            description="Number of atoms in the structure, equal to the size of the Point Domain",
        )
        index_of_first = tree.outputs.integer(
            "Index of First",
            description="`Index` of first atom in the entire structure, which will always be 0",
        )
        index_of_last = tree.outputs.integer(
            "Index of Last",
            description="Index of last atom in the entire structure. Equal to `Atom Count` - 1",
        )
        residue_factor = tree.outputs.float("Residue Factor", subtype="FACTOR")
        residue_index = tree.outputs.integer(
            "Residue Index",
            description="A unique `Group ID` that increases whenever `Sub Group ID` or `Group ID` change",
        )
        residue_count = tree.outputs.integer("Residue Count")

        group = SubGroupInfo(sub_group_id=UniqueResidueID())
        (
            IndexToFactor(index=group.o.group_id, size=group.o.sub_group_total)
            >> residue_factor
        )
        group_1 = GroupInfo()
        index_1 = g.Index()
        IndexToFactor(index=index_1, size=group_1.o.size) >> atom_factor

        index_1 >> index
        group_1 >> atom_count
        group_1.o.index_of_first >> index_of_first
        group_1.o.index_of_last >> index_of_last
        group.o.group_id >> residue_index
        group.o.sub_group_total >> residue_count


ASSET = StructureParameter

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
