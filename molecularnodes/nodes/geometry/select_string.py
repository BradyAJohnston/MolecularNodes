# Node-group asset "Select String" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    CustomGeometryGroup,
    PackageLibrary,
    SocketAccessor,
    StringSocket,
)
from nodebpy.types import InputString
from .chain_id import ChainID
from .residue_id import ResidueID
from .residue_name import ResidueName


class StringListIndex(CustomGeometryGroup):
    _name = "String List Index"
    _color_tag = "CONVERTER"
    _tree_properties = {
        "description": "Position of a string in a comma-separated list. Each entry may hold several aliases separated by '/', any of which matches. -1 when not found or the item is empty."
    }

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        list = tree.inputs.string(
            "List", "", description="Comma-separated entries to search"
        )
        item = tree.inputs.string("Item", "", description="Entry to look for")
        index = tree.outputs.integer(
            "Index",
            description="Position of the item in the list, or -1 when it is not found",
        )

        trim_string = list.split(",").trim()
        repeat_zone = g.RepeatZone(trim_string.list_length())
        index_1 = repeat_zone.items.integer("Index", -1)
        with g.Frame("Match entry"):
            _string = g.String(
                string="Wrapping both the entry and the item in '/' makes the match exact while still letting an entry list aliases like MET/MSE. The item is trimmed first, which also keeps the group input off the Join Strings multi-input, whose link order the asset build does not preserve."
            )
            join_strings = g.JoinStrings(
                (
                    g.JoinStrings(
                        (g.String(string="/"), trim_string[repeat_zone.iteration])
                    ),
                    g.String(string="/"),
                )
            )
            join_strings_1 = g.JoinStrings(
                (
                    g.JoinStrings((g.String(string="/"), item.trim())),
                    g.String(string="/"),
                )
            )
            switch = join_strings.o.string.contains(join_strings_1).switch.integer(
                index_1.current, repeat_zone.iteration
            )
        switch >> index_1.next
        (
            g.Compare.integer.equal(item.length(), 0).o.result.switch.integer(
                index_1.result, -1
            )
            >> index
        )


class SelectString(AssetGeometryGroup):
    """
    Select atoms with a comma-separated list of chain IDs, residue IDs or ranges, and residue names. 'A:10-20' restricts a chain to residues, 'A:' selects a chain whose ID looks like a number and ':A' a residue name that matches a chain.

    Parameters
    ----------
    selection : InputString
        Comma-separated terms, any of which selects an atom: a chain ID (A), a residue ID (42) or range (10-20), a residue name (LYS), or a chain with residues after a colon (A:10-20, A:LYS). Use 'A:' for a chain named like a number and ':A' for a residue name that is also a chain ID.
    chain_ids : InputString
        Comma-separated chain IDs in the order they are numbered in the `chain_id` attribute. Filled in from the structure when the node is added to an entity's tree.
    residue_names : InputString
        Comma-separated residue names in the order they are numbered in the `res_name` attribute. Names sharing a number are separated by '/'.

    Inputs
    ------
    i.selection : StringSocket
        Comma-separated terms, any of which selects an atom: a chain ID (A), a residue ID (42) or range (10-20), a residue name (LYS), or a chain with residues after a colon (A:10-20, A:LYS). Use 'A:' for a chain named like a number and ':A' for a residue name that is also a chain ID.
    i.chain_ids : StringSocket
        Comma-separated chain IDs in the order they are numbered in the `chain_id` attribute. Filled in from the structure when the node is added to an entity's tree.
    i.residue_names : StringSocket
        Comma-separated residue names in the order they are numbered in the `res_name` attribute. Names sharing a number are separated by '/'.

    Outputs
    -------
    o.selection : BooleanSocket
        Selection
    o.inverted : BooleanSocket
        Inverted
    """

    _name = "Select String"
    _asset_name = "Select String"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"
    _tree_properties = {
        "description": "Select atoms with a comma-separated list of chain IDs, residue IDs or ranges, and residue names. 'A:10-20' restricts a chain to residues, 'A:' selects a chain whose ID looks like a number and ':A' a residue name that matches a chain."
    }

    class _Inputs(SocketAccessor):
        selection: StringSocket
        """Comma-separated terms, any of which selects an atom: a chain ID (A), a residue ID (42) or range (10-20), a residue name (LYS), or a chain with residues after a colon (A:10-20, A:LYS). Use 'A:' for a chain named like a number and ':A' for a residue name that is also a chain ID."""
        chain_ids: StringSocket
        """Comma-separated chain IDs in the order they are numbered in the `chain_id` attribute. Filled in from the structure when the node is added to an entity's tree."""
        residue_names: StringSocket
        """Comma-separated residue names in the order they are numbered in the `res_name` attribute. Names sharing a number are separated by '/'."""

    class _Outputs(SocketAccessor):
        selection: BooleanSocket
        """Selection"""
        inverted: BooleanSocket
        """Inverted"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        selection: InputString = "",
        chain_ids: InputString = "",
        residue_names: InputString = "ALA,ARG,ASN,ASP/ASH,CYS/CYM/CYX,GLU/GLH,GLN,GLY,HIS/HID/HIE/HIP/HYP,ILE,LEU,LYS/LYN,MET/MSE,PHE,PRO,SER/SNC,THR,TRP,TYR,VAL,,,,,,,,,,,DA,DC,DG,DT/PST,,,,,,,A,C,G,U/T",
    ):
        super().__init__(
            **{
                "Selection": selection,
                "Chain IDs": chain_ids,
                "Residue Names": residue_names,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        selection = tree.inputs.string(
            "Selection",
            "",
            description="Comma-separated terms, any of which selects an atom: a chain ID (A), a residue ID (42) or range (10-20), a residue name (LYS), or a chain with residues after a colon (A:10-20, A:LYS). Use 'A:' for a chain named like a number and ':A' for a residue name that is also a chain ID.",
            optional_label=True,
        )
        with tree.inputs.panel("Context", default_closed=True):
            chain_ids = tree.inputs.string(
                "Chain IDs",
                "",
                description="Comma-separated chain IDs in the order they are numbered in the `chain_id` attribute. Filled in from the structure when the node is added to an entity's tree.",
            )
            residue_names = tree.inputs.string(
                "Residue Names",
                "ALA,ARG,ASN,ASP/ASH,CYS/CYM/CYX,GLU/GLH,GLN,GLY,HIS/HID/HIE/HIP/HYP,ILE,LEU,LYS/LYN,MET/MSE,PHE,PRO,SER/SNC,THR,TRP,TYR,VAL,,,,,,,,,,,DA,DC,DG,DT/PST,,,,,,,A,C,G,U/T",
                description="Comma-separated residue names in the order they are numbered in the `res_name` attribute. Names sharing a number are separated by '/'.",
            )
        selection_1 = tree.outputs.boolean("Selection")
        inverted = tree.outputs.boolean("Inverted")

        trim_string = selection.split(",").trim()
        repeat_zone = g.RepeatZone(trim_string.list_length())
        selection_2 = repeat_zone.items.boolean("Selection")
        with g.Frame("Split term"):
            _string = g.String(
                string="A term is 'chain:residue'. Without a colon the whole term is tried as a residue ID or range first, then as a chain ID, then as a residue name."
            )
            get_list_item = trim_string[repeat_zone.iteration]
            match_string = get_list_item.contains(":")
            switch = match_string.switch.string(
                g.JoinStrings((g.String(string=":"), get_list_item)), get_list_item
            )
            split_string = switch.split(":")
            trim_string_1 = split_string[0].trim()
            trim_string_2 = split_string[1].trim()
        with g.Frame("Chain"):
            string_list_index = StringListIndex(
                List=chain_ids,
                Item=match_string.switch.string(trim_string_2, trim_string_1),
            )
            compare = string_list_index >= 0
            boolean_math = compare & g.Compare.integer.equal(
                ChainID(), string_list_index
            )
        with g.Frame("Residue ID"):
            residue_id = ResidueID()
            _string_1 = g.String(
                string="A '-' after the first character marks a range, so a negative residue ID like -5 still parses as a single number."
            )
            split_string_1 = trim_string_2.split("-")
            match_string_1 = trim_string_2.slice(1, trim_string_2.length()).contains(
                "-"
            )
            boolean_math_1 = (residue_id >= split_string_1[0].to_integer()) & (
                residue_id <= split_string_1[1].to_integer()
            )
            string_to_value = trim_string_2.to_integer()
            compare_1 = g.Compare.integer.equal(residue_id, string_to_value)
            compare_2 = g.Compare.string.equal(
                string_to_value.to_string(), trim_string_2
            )
        with g.Frame("Residue name"):
            string_list_index_1 = StringListIndex(
                List=residue_names, Item=trim_string_2
            )
            boolean_math_2 = (string_list_index_1 >= 0) & g.Compare.integer.equal(
                ResidueName(), string_list_index_1
            )
        with g.Frame("Combine"):
            switch_1 = match_string_1.switch.boolean(
                compare_2.o.result.switch.boolean(boolean_math_2, compare_1),
                boolean_math_1,
            )
            boolean_math_3 = (
                g.Compare.integer.equal(trim_string_1.length(), 0).o.result
                | boolean_math
            ) & (g.Compare.integer.equal(trim_string_2.length(), 0).o.result | switch_1)
            switch_2 = match_string.switch.boolean(
                (compare & ~(match_string_1 | compare_2)).switch.boolean(
                    switch_1, boolean_math
                ),
                boolean_math_3,
            )
            boolean_math_4 = selection_2.current | switch_2
        boolean_math_4 >> selection_2.next
        ~selection_2.result >> inverted

        selection_2.result >> selection_1


ASSET = SelectString

ASSET_METADATA = {
    "description": "Select atoms with a comma-separated list of chain IDs, residue IDs or ranges, and residue names",
    "catalog_id": "bd1f205b-fea5-4700-b2c2-754f3321e969",
}
