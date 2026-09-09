# Node-group asset 'Select Res ID String' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    PackageLibrary,
    SocketAccessor,
    StringSocket,
)
from nodebpy.types import InputString
from .residue_id import ResidueID


class SelectResIDString(AssetGeometryGroup):
    """
    Select Res ID String

    Parameters
    ----------
    string : InputString
        String

    Inputs
    ------
    i.string : StringSocket
        String

    Outputs
    -------
    o.selection : BooleanSocket
        Selection
    """

    _name = "Select Res ID String"
    _asset_name = "Select Res ID String"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        string: StringSocket
        """String"""

    class _Outputs(SocketAccessor):
        selection: BooleanSocket
        """Selection"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        string: InputString = "",
    ):
        super().__init__(**{"String": string})

    def _build_group(self, tree):
        string = tree.inputs.string("String", "", optional_label=True)
        selection = tree.outputs.boolean("Selection")

        closure_zone = g.ClosureZone()
        string_1 = closure_zone.inputs.string("String")
        selection_1 = closure_zone.outputs.boolean("Selection")
        closure_zone_1 = g.ClosureZone()
        string_2 = closure_zone_1.inputs.string("String")
        selection_2 = closure_zone_1.outputs.boolean("Selection")
        group = ResidueID()
        string_3 = g.String(string="-")
        (
            g.Compare.integer.equal(
                g.StringToValue.integer(string_1).o.value, ResidueID()
            )
            >> selection_1
        )
        trim_string = g.TrimString(string=g.SplitString(string=string, separator=","))
        repeat_zone = g.RepeatZone(trim_string.o.string.list_length())
        boolean = repeat_zone.items.boolean("Boolean")
        get_list_item = trim_string.o.string[repeat_zone.iteration]
        trim_string_1 = g.TrimString(
            string=g.SplitString(string=string_2, separator=string_3)
        )
        boolean_math = (
            group >= g.StringToValue.integer(trim_string_1.o.string[0]).o.value
        ) & (group <= g.StringToValue.integer(trim_string_1.o.string[1]).o.value)
        boolean_math >> selection_2
        switch = get_list_item.contains(string_3).switch.closure(
            closure_zone.closure, closure_zone_1.closure
        )
        evaluate_closure = g.EvaluateClosure(switch)
        evaluate_closure.inputs.string("String", get_list_item)
        selection_3 = evaluate_closure.outputs.boolean("Selection")
        (boolean.current | selection_3) >> boolean.next

        boolean.result >> selection


ASSET = SelectResIDString

ASSET_METADATA = {
    "catalog_id": "bd1f205b-fea5-4700-b2c2-754f3321e969",
}
