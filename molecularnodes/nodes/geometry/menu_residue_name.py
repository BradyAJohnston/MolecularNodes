# Node-group asset "Menu Residue Name" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    IntegerSocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputMenu
from .residue_name import ResidueName


class MenuResidueName(AssetGeometryGroup):
    """
    Menu Residue Name

    Parameters
    ----------
    residue_name : InputMenu | Literal["UNK", "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "SNC", "MSE", "ASH", "CYM", "CYX", "GLH", "HID", "HIE", "HIP", "HYP", "LYN", "DA", "DC", "DG", "DT", "PST", "rA", "rC", "rG", "rU"]
        Residue Name

    Inputs
    ------
    i.residue_name : MenuSocket
        Residue Name

    Outputs
    -------
    o.res_name : IntegerSocket
        res_name
    o.selection : BooleanSocket
        Selection
    """

    _name = "Menu Residue Name"
    _asset_name = "Menu Residue Name"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        residue_name: MenuSocket
        """Residue Name"""

    class _Outputs(SocketAccessor):
        res_name: IntegerSocket
        selection: BooleanSocket
        """Selection"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        residue_name: InputMenu
        | Literal[
            "UNK",
            "ALA",
            "ARG",
            "ASN",
            "ASP",
            "CYS",
            "GLU",
            "GLN",
            "GLY",
            "HIS",
            "ILE",
            "LEU",
            "LYS",
            "MET",
            "PHE",
            "PRO",
            "SER",
            "THR",
            "TRP",
            "TYR",
            "VAL",
            "SNC",
            "MSE",
            "ASH",
            "CYM",
            "CYX",
            "GLH",
            "HID",
            "HIE",
            "HIP",
            "HYP",
            "LYN",
            "DA",
            "DC",
            "DG",
            "DT",
            "PST",
            "rA",
            "rC",
            "rG",
            "rU",
        ] = "ALA",
    ):
        super().__init__(**{"Residue Name": residue_name})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        residue_name = tree.inputs.menu("Residue Name", optional_label=True)
        res_name = tree.outputs.integer("res_name")
        selection = tree.outputs.boolean("Selection")

        menu_switch = g.MenuSwitch.integer(
            residue_name,
            {
                "UNK": -1,
                "ALA": 0,
                "ARG": 1,
                "ASN": 2,
                "ASP": 3,
                "CYS": 4,
                "GLU": 5,
                "GLN": 6,
                "GLY": 7,
                "HIS": 8,
                "ILE": 9,
                "LEU": 10,
                "LYS": 11,
                "MET": 12,
                "PHE": 13,
                "PRO": 14,
                "SER": 15,
                "THR": 16,
                "TRP": 17,
                "TYR": 18,
                "VAL": 19,
                "SNC": 15,
                "MSE": 12,
                "ASH": 3,
                "CYM": 4,
                "CYX": 4,
                "GLH": 5,
                "HID": 8,
                "HIE": 8,
                "HIP": 8,
                "HYP": 8,
                "LYN": 11,
                "DA": 30,
                "DC": 31,
                "DG": 32,
                "DT": 33,
                "PST": 33,
                "rA": 40,
                "rC": 41,
                "rG": 42,
                "rU": 43,
            },
        )
        g.Compare.integer.equal(menu_switch.o.output, ResidueName()) >> selection

        menu_switch >> res_name

        residue_name.default_value = "ALA"


ASSET = MenuResidueName

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
