# Node-group asset 'Menu Residue Mask' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    IntegerSocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputMenu
from .group_pick_vector import GroupPickVector
from .menu_atom_name import MenuAtomName
from .ures_id import UResID


class MenuResidueMask(AssetGeometryGroup):
    """
    Menu Residue Mask

    Parameters
    ----------
    atom_name : InputMenu | Literal["N", "CA", "C", "O", "CB", "CG", "CG1", "CG2", "OG", "OG1", "SG", "CD", "CD1", "CD2", "ND1", "ND2", "OD1", "OD2", "SD", "CE", "CE1", "CE2", "CE3", "NE", "NE1", "NE2", "OE1", "OE2", "CH2", "NH1", "NH2", "OH", "CZ", "CZ2", "CZ3", "NZ", "OXT", "P", "O1P", "OP1", "OP2", "O2P", "O5", "C5", "C4", "O4", "C3", "O3", "C2", "O2", "C1", "N1", "N9", "N3", "C8", "N7", "C5", "C6", "N6", "C2", "C4", "O6", "N2", "N4", "O2", "O4", "C7"]
        Atom Name

    Inputs
    ------
    i.atom_name : MenuSocket
        Atom Name

    Outputs
    -------
    o.is_valid : BooleanSocket
        Group contains only one occurrance of the selected atom. None or more than one returns False
    o.index : IntegerSocket
        Index for the group's atom with specified name, returns -1 if not valid
    o.position : VectorSocket
        Position of the picked point in the group, returns (0, 0, 0) if not valid
    """

    _name = "Menu Residue Mask"
    _asset_name = "Menu Residue Mask"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.menu_residue_mask"}

    class _Inputs(SocketAccessor):
        atom_name: MenuSocket
        """Atom Name"""

    class _Outputs(SocketAccessor):
        is_valid: BooleanSocket
        """Group contains only one occurrance of the selected atom. None or more than one returns False"""
        index: IntegerSocket
        """Index for the group's atom with specified name, returns -1 if not valid"""
        position: VectorSocket
        """Position of the picked point in the group, returns (0, 0, 0) if not valid"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        atom_name: InputMenu
        | Literal[
            "N",
            "CA",
            "C",
            "O",
            "CB",
            "CG",
            "CG1",
            "CG2",
            "OG",
            "OG1",
            "SG",
            "CD",
            "CD1",
            "CD2",
            "ND1",
            "ND2",
            "OD1",
            "OD2",
            "SD",
            "CE",
            "CE1",
            "CE2",
            "CE3",
            "NE",
            "NE1",
            "NE2",
            "OE1",
            "OE2",
            "CH2",
            "NH1",
            "NH2",
            "OH",
            "CZ",
            "CZ2",
            "CZ3",
            "NZ",
            "OXT",
            "P",
            "O1P",
            "OP1",
            "OP2",
            "O2P",
            "O5",
            "C5",
            "C4",
            "O4",
            "C3",
            "O3",
            "C2",
            "O2",
            "C1",
            "N1",
            "N9",
            "N3",
            "C8",
            "N7",
            "C5",
            "C6",
            "N6",
            "C2",
            "C4",
            "O6",
            "N2",
            "N4",
            "O2",
            "O4",
            "C7",
        ] = "N",
    ):
        super().__init__(**{"Atom Name": atom_name})

    def _build_group(self, tree):
        atom_name = tree.inputs.menu("Atom Name", optional_label=True)
        is_valid = tree.outputs.boolean(
            "Is Valid",
            description="Group contains only one occurrance of the selected atom. None or more than one returns False",
        )
        index = tree.outputs.integer(
            "Index",
            description="Index for the group's atom with specified name, returns -1 if not valid",
        )
        position = tree.outputs.vector(
            "Position",
            description="Position of the picked point in the group, returns (0, 0, 0) if not valid",
        )

        group = GroupPickVector(
            pick=MenuAtomName(atom_name=atom_name).o.selection,
            group_id=UResID().o.ures_id,
        )

        group >> is_valid
        group.o.index >> index
        group.o.vector >> position

        atom_name.default_value = "N"


ASSET = MenuResidueMask

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
