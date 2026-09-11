# Node-group asset "Find Bonded Atom" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
    VectorSocket,
)
from nodebpy.types import InputInteger, InputMenu
from ._shared.find_connected import FindConnected
from .atom_name import AtomName
from .menu_atom_name import MenuAtomName


class FindBondedAtom(AssetGeometryGroup):
    """
    Find Bonded Atom

    Parameters
    ----------
    index : InputInteger
        Index
    method : InputMenu | Literal["Any", "Exact"]
        Method
    atom_name : InputMenu | Literal["N", "CA", "C", "O", "CB", "CG", "CG1", "CG2", "OG", "OG1", "SG", "CD", "CD1", "CD2", "ND1", "ND2", "OD1", "OD2", "SD", "CE", "CE1", "CE2", "CE3", "NE", "NE1", "NE2", "OE1", "OE2", "CH2", "NH1", "NH2", "OH", "CZ", "CZ2", "CZ3", "NZ", "OXT", "P", "O1P", "OP1", "OP2", "O2P", "O5", "C5", "C4", "O4", "C3", "O3", "C2", "O2", "C1", "N1", "N9", "N3", "C8", "N7", "C5", "C6", "N6", "C2", "C4", "O6", "N2", "N4", "O2", "O4", "C7"]
        Atom Name
    distance : InputInteger
        Distance

    Inputs
    ------
    i.index : IntegerSocket
        Index
    i.method : MenuSocket
        Method
    i.atom_name : MenuSocket
        Atom Name
    i.distance : IntegerSocket
        Distance

    Outputs
    -------
    o.is_valid : BooleanSocket
        Is Valid
    o.index : IntegerSocket
        Index
    o.position : VectorSocket
        Position
    """

    _name = "Find Bonded Atom"
    _asset_name = "Find Bonded Atom"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        index: IntegerSocket
        """Index"""
        method: MenuSocket
        """Method"""
        atom_name: MenuSocket
        """Atom Name"""
        distance: IntegerSocket
        """Distance"""

    class _Outputs(SocketAccessor):
        is_valid: BooleanSocket
        """Is Valid"""
        index: IntegerSocket
        """Index"""
        position: VectorSocket
        """Position"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        index: InputInteger = 0,
        method: InputMenu | Literal["Any", "Exact"] = "Exact",
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
        distance: InputInteger = 2,
    ):
        super().__init__(
            **{
                "Index": index,
                "Method": method,
                "Atom Name": atom_name,
                "Distance": distance,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        index = tree.inputs.integer("Index", 0, min_value=0, default_input="INDEX")
        method = tree.inputs.menu("Method", optional_label=True)
        atom_name = tree.inputs.menu("Atom Name", optional_label=True)
        distance = tree.inputs.integer("Distance", 2, min_value=1, max_value=3)
        is_valid = tree.outputs.boolean("Is Valid")
        index_1 = tree.outputs.integer("Index")
        position = tree.outputs.vector("Position")

        group = FindConnected(
            value=AtomName(),
            match=MenuAtomName(atom_name=atom_name).o[0],
            distance=distance,
            method=method,
        )
        evaluate_at_index = group.o.index.point.at(index)
        g.Position().o.position.point.at(evaluate_at_index) >> position

        group >> is_valid
        evaluate_at_index >> index_1

        method.default_value = "Exact"
        atom_name.default_value = "N"


ASSET = FindBondedAtom

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
