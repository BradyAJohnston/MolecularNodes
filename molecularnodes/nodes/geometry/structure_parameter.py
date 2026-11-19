# Node-group asset "Structure Parameter" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
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
from .ures_id import UResID


class StructureParameter(AssetGeometryGroup):
    """
    Structure Parameter

    Outputs
    -------
    o.factor : FloatSocket
        An atom's relative position in the structure, 0 for the first atom and 1 for the last
    o.atom_index : IntegerSocket
        Index of the atom within the structure, equal to the `Index` input node
    o.atom_count : IntegerSocket
        Number of atoms in the structure
    o.residue_factor : FloatSocket
        A residue's relative position in the structure, 0 for the first residue and 1 for the last
    o.residue_index : IntegerSocket
        Index of the residue within the structure, starting at 0 for the first residue
    o.residue_count : IntegerSocket
        Number of residues in the structure
    """

    _name = "Structure Parameter"
    _asset_name = "Structure Parameter"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        factor: FloatSocket
        """An atom's relative position in the structure, 0 for the first atom and 1 for the last"""
        atom_index: IntegerSocket
        """Index of the atom within the structure, equal to the `Index` input node"""
        atom_count: IntegerSocket
        """Number of atoms in the structure"""
        residue_factor: FloatSocket
        """A residue's relative position in the structure, 0 for the first residue and 1 for the last"""
        residue_index: IntegerSocket
        """Index of the residue within the structure, starting at 0 for the first residue"""
        residue_count: IntegerSocket
        """Number of residues in the structure"""

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
            description="An atom's relative position in the structure, 0 for the first atom and 1 for the last",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        atom_index = tree.outputs.integer(
            "Atom Index",
            description="Index of the atom within the structure, equal to the `Index` input node",
        )
        atom_count = tree.outputs.integer(
            "Atom Count", description="Number of atoms in the structure"
        )
        residue_factor = tree.outputs.float(
            "Residue Factor",
            description="A residue's relative position in the structure, 0 for the first residue and 1 for the last",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        residue_index = tree.outputs.integer(
            "Residue Index",
            description="Index of the residue within the structure, starting at 0 for the first residue",
        )
        residue_count = tree.outputs.integer(
            "Residue Count", description="Number of residues in the structure"
        )

        group_info = GroupInfo()
        index = g.Index()
        IndexToFactor(index=index, size=group_info.o.size) >> factor
        sub_group_info = SubGroupInfo(sub_group_id=UResID().o.ures_id)
        (
            IndexToFactor(
                index=sub_group_info.o.group_id, size=sub_group_info.o.sub_group_total
            )
            >> residue_factor
        )

        index >> atom_index
        group_info >> atom_count
        sub_group_info.o.group_id >> residue_index
        sub_group_info.o.sub_group_total >> residue_count


ASSET = StructureParameter

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
