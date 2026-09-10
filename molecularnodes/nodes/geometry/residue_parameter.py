# Node-group asset "Residue Parameter" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from ._shared.index_to_factor import IndexToFactor
from .atom_name import AtomName
from .chain_id import ChainID
from .residue_id import ResidueID
from .sub_group_info import SubGroupInfo


class ResidueParameter(AssetGeometryGroup):
    """
    Residue Parameter

    Outputs
    -------
    o.factor : FloatSocket
        An atom's relative position in a residue, with the first atom being 0 and the last atom being 1
    o.atom_count : IntegerSocket
        Number of atoms in a residue
    o.atom_index : IntegerSocket
        Index of an atom in a residue when counting from 0
    o.first_atom_name : IntegerSocket
        the atom_name for the first atom in a residue
    o.last_atom_name : IntegerSocket
        The atom_name for the last atom in a residue
    o.index_of_first : IntegerSocket
        Index (in the whole structure) for the first atom in a residue
    o.index_of_last : IntegerSocket
        Index (in the whole structure) for the last atom in a residue
    """

    _name = "Residue Parameter"
    _asset_name = "Residue Parameter"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.residue_parameter"}

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        factor: FloatSocket
        """An atom's relative position in a residue, with the first atom being 0 and the last atom being 1"""
        atom_count: IntegerSocket
        """Number of atoms in a residue"""
        atom_index: IntegerSocket
        """Index of an atom in a residue when counting from 0"""
        first_atom_name: IntegerSocket
        """the atom_name for the first atom in a residue"""
        last_atom_name: IntegerSocket
        """The atom_name for the last atom in a residue"""
        index_of_first: IntegerSocket
        """Index (in the whole structure) for the first atom in a residue"""
        index_of_last: IntegerSocket
        """Index (in the whole structure) for the last atom in a residue"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self):
        super().__init__()

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        factor = tree.outputs.float(
            "Factor",
            description="An atom's relative position in a residue, with the first atom being 0 and the last atom being 1",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        atom_count = tree.outputs.integer(
            "Atom Count", description="Number of  atoms in a residue"
        )
        atom_index = tree.outputs.integer(
            "Atom Index",
            description="Index of an atom in a residue when counting from 0",
        )
        first_atom_name = tree.outputs.integer(
            "First atom_name",
            description="the atom_name for the first atom in a residue",
        )
        last_atom_name = tree.outputs.integer(
            "Last atom_name", description="The atom_name for the last atom in a residue"
        )
        index_of_first = tree.outputs.integer(
            "Index of First",
            description="Index (in the whole structure) for the first atom in a  residue",
        )
        index_of_last = tree.outputs.integer(
            "Index of Last",
            description="Index (in the whole structure) for the last atom in a  residue",
        )

        group = SubGroupInfo(sub_group_id=ResidueID(), group_id=ChainID())
        IndexToFactor(index=group.o.index_in_group_id, size=group.o.size) >> factor
        AtomName(index=group.o.index_of_first) >> first_atom_name
        AtomName(index=group.o.index_of_last) >> last_atom_name

        group >> atom_count
        group.o.index_in_group_id >> atom_index
        group.o.index_of_first >> index_of_first
        group.o.index_of_last >> index_of_last


ASSET = ResidueParameter

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
