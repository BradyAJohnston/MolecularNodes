# Node-group asset 'Menu Atom Name' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    IntegerSocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
    StringSocket,
)
from nodebpy.types import InputMenu
from .atom_name import AtomName


class MenuAtomName(AssetGeometryGroup):
    """
    Menu Atom Name

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
    o.socket_1 : IntegerSocket
        atom_name
    o.selection : BooleanSocket
        Selection
    o.socket_3 : StringSocket
        atom_name
    """

    _name = "Menu Atom Name"
    _asset_name = "Menu Atom Name"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        atom_name: MenuSocket
        """Atom Name"""

    class _Outputs(SocketAccessor):
        socket_1: IntegerSocket
        """atom_name"""
        selection: BooleanSocket
        """Selection"""
        socket_3: StringSocket
        """atom_name"""

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
        atom_name_1 = tree.outputs.integer("atom_name")
        selection = tree.outputs.boolean("Selection")
        atom_name_2 = tree.outputs.string("atom_name")

        _index_switch = g.IndexSwitch.integer(
            items=(
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8,
                9,
                10,
                11,
                12,
                13,
                14,
                15,
                16,
                17,
                18,
                19,
                20,
                21,
                23,
                24,
                25,
                26,
                27,
                28,
                29,
                30,
                31,
                32,
                33,
                34,
                35,
                36,
                37,
                38,
                50,
                51,
                51,
                52,
                52,
                53,
                54,
                55,
                56,
                57,
                58,
                59,
                60,
                61,
                62,
                63,
                64,
                65,
                66,
                67,
                68,
                69,
                70,
                71,
                72,
                73,
                74,
                75,
                76,
                77,
            )
        )
        menu_switch = g.MenuSwitch.integer(
            atom_name,
            {
                "N": 1,
                "CA": 2,
                "C": 3,
                "O": 4,
                "CB": 5,
                "CG": 6,
                "CG1": 7,
                "CG2": 8,
                "OG": 9,
                "OG1": 10,
                "SG": 11,
                "CD": 12,
                "CD1": 13,
                "CD2": 14,
                "ND1": 15,
                "ND2": 16,
                "OD1": 17,
                "OD2": 18,
                "SD": 19,
                "CE": 20,
                "CE1": 21,
                "CE2": 23,
                "CE3": 24,
                "NE": 25,
                "NE1": 26,
                "NE2": 27,
                "OE1": 28,
                "OE2": 29,
                "CH2": 30,
                "NH1": 31,
                "NH2": 32,
                "OH": 33,
                "CZ": 34,
                "CZ2": 35,
                "CZ3": 36,
                "NZ": 37,
                "OXT": 38,
                "P": 50,
                "O1P": 51,
                "OP1": 51,
                "OP2": 52,
                "O2P": 52,
                "O5'": 53,
                "C5'": 54,
                "C4'": 55,
                "O4'": 56,
                "C3'": 57,
                "O3'": 58,
                "C2'": 59,
                "O2'": 60,
                "C1'": 61,
                "N1": 62,
                "N9": 63,
                "N3": 64,
                "C8": 65,
                "N7": 66,
                "C5": 67,
                "C6": 68,
                "N6": 69,
                "C2": 70,
                "C4": 71,
                "O6": 72,
                "N2": 73,
                "N4": 74,
                "O2": 75,
                "O4": 76,
                "C7": 77,
            },
        )
        g.Compare.integer.equal(menu_switch.o.output, AtomName()) >> selection
        (
            g.IndexSwitch.string(
                menu_switch.o.output,
                (
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
                    "O5'",
                    "C5'",
                    "C4'",
                    "O4'",
                    "C3'",
                    "O3'",
                    "C2'",
                    "O2'",
                    "C1'",
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
                ),
            )
            >> atom_name_2
        )
        _menu_switch_1 = g.MenuSwitch.string(
            items={
                "N": "N",
                "CA": "CA",
                "C": "C",
                "O": "O",
                "CB": "CB",
                "CG": "CG",
                "CG1": "CG1",
                "CG2": "CG2",
                "OG": "OG",
                "OG1": "OG1",
                "SG": "SG",
                "CD": "CD",
                "CD1": "CD1",
                "CD2": "CD2",
                "ND1": "ND1",
                "ND2": "ND2",
                "OD1": "OD1",
                "OD2": "OD2",
                "SD": "SD",
                "CE": "CE",
                "CE1": "CE1",
                "CE2": "CE2",
                "CE3": "CE3",
                "NE": "NE",
                "NE1": "NE1",
                "NE2": "NE2",
                "OE1": "OE1",
                "OE2": "OE2",
                "CH2": "CH2",
                "NH1": "NH1",
                "NH2": "NH2",
                "OH": "OH",
                "CZ": "CZ",
                "CZ2": "CZ2",
                "CZ3": "CZ3",
                "NZ": "NZ",
                "OXT": "OXT",
                "P": "P",
                "O1P": "O1P",
                "OP1": "OP1",
                "OP2": "OP2",
                "O2P": "O2P",
                "O5'": "O5'",
                "C5'": "C5'",
                "C4'": "C4'",
                "O4'": "O4'",
                "C3'": "C3'",
                "O3'": "O3'",
                "C2'": "C2'",
                "O2'": "O2'",
                "C1'": "C1'",
                "N1": "N1",
                "N9": "N9",
                "N3": "N3",
                "C8": "C8",
                "N7": "N7",
                "C5": "C5",
                "C6": "C6",
                "N6": "N6",
                "C2": "C2",
                "C4": "C4",
                "O6": "O6",
                "N2": "N2",
                "N4": "N4",
                "O2": "O2",
                "O4": "O4",
                "C7": "C7",
            }
        )

        menu_switch >> atom_name_1

        atom_name.default_value = "N"


ASSET = MenuAtomName

ASSET_METADATA = {
    "catalog_id": "bd1f205b-fea5-4700-b2c2-754f3321e969",
}
