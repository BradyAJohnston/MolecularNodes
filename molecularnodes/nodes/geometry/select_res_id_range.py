# Node-group asset 'Select Res ID Range' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputInteger
from .boolean_andor import BooleanAndOr
from .residue_id import ResidueID


class SelectResIDRange(AssetGeometryGroup):
    """
    Select Res ID Range

    Parameters
    ----------
    and_ : InputBoolean
        The resulting selection must overlap with this input selection
    or_ : InputBoolean
        The resulting selection can be calculated from this node or be from this input selection
    min : InputInteger
        Minimum of a `res_id` range selection
    max : InputInteger
        Maximum of a `res_id` range selection

    Inputs
    ------
    i.and_ : BooleanSocket
        The resulting selection must overlap with this input selection
    i.or_ : BooleanSocket
        The resulting selection can be calculated from this node or be from this input selection
    i.min : IntegerSocket
        Minimum of a `res_id` range selection
    i.max : IntegerSocket
        Maximum of a `res_id` range selection

    Outputs
    -------
    o.selection : BooleanSocket
        The calculated selection
    o.inverted : BooleanSocket
        The inverse of the calculated selection
    """

    _name = "Select Res ID Range"
    _asset_name = "Select Res ID Range"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.select_res_id_range"}

    class _Inputs(SocketAccessor):
        and_: BooleanSocket
        """The resulting selection must overlap with this input selection"""
        or_: BooleanSocket
        """The resulting selection can be calculated from this node or be from this input selection"""
        min: IntegerSocket
        """Minimum of a `res_id` range selection"""
        max: IntegerSocket
        """Maximum of a `res_id` range selection"""

    class _Outputs(SocketAccessor):
        selection: BooleanSocket
        """The calculated selection"""
        inverted: BooleanSocket
        """The inverse of the calculated selection"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        and_: InputBoolean = True,
        or_: InputBoolean = False,
        min: InputInteger = 10,
        max: InputInteger = 100,
    ):
        super().__init__(**{"And": and_, "Or": or_, "Min": min, "Max": max})

    def _build_group(self, tree):
        and_ = tree.inputs.boolean(
            "And",
            True,
            description="The resulting selection must overlap with this input selection",
            hide_value=True,
        )
        or_ = tree.inputs.boolean(
            "Or",
            False,
            description="The resulting selection can be calculated from this node or be from this input selection",
            hide_value=True,
        )
        with tree.inputs.panel("Res ID"):
            min = tree.inputs.integer(
                "Min",
                10,
                description="Minimum of a `res_id` range selection",
                min_value=0,
            )
            max = tree.inputs.integer(
                "Max",
                100,
                description="Maximum of a `res_id` range selection",
                min_value=1,
            )
        selection = tree.outputs.boolean(
            "Selection", description="The calculated selection"
        )
        inverted = tree.outputs.boolean(
            "Inverted", description="The inverse of the calculated selection"
        )

        group = ResidueID()
        group_1 = BooleanAndOr(
            and_=and_, or_=or_, boolean=(group >= min) & (group <= max)
        )

        group_1 >> selection
        group_1.o.inverted >> inverted


ASSET = SelectResIDRange

ASSET_METADATA = {
    "catalog_id": "bd1f205b-fea5-4700-b2c2-754f3321e969",
}
