# Node-group asset "Select Nucleic Type" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean
from .boolean_andor import BooleanAndOr
from .residue_name import ResidueName


class SelectNucleicType(AssetGeometryGroup):
    """
    Select Nucleic Type

    Parameters
    ----------
    and_ : InputBoolean
        The resulting selection must overlap with this input selection
    or_ : InputBoolean
        The resulting selection can be calculated from this node or be from this input selection

    Inputs
    ------
    i.and_ : BooleanSocket
        The resulting selection must overlap with this input selection
    i.or_ : BooleanSocket
        The resulting selection can be calculated from this node or be from this input selection

    Outputs
    -------
    o.is_purine : BooleanSocket
        is_purine
    o.is_pyrimidine : BooleanSocket
        is_pyrimidine
    """

    _name = "Select Nucleic Type"
    _asset_name = "Select Nucleic Type"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.select_nucleic_type"}

    class _Inputs(SocketAccessor):
        and_: BooleanSocket
        """The resulting selection must overlap with this input selection"""
        or_: BooleanSocket
        """The resulting selection can be calculated from this node or be from this input selection"""

    class _Outputs(SocketAccessor):
        is_purine: BooleanSocket
        is_pyrimidine: BooleanSocket

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        and_: InputBoolean = True,
        or_: InputBoolean = False,
    ):
        super().__init__(**{"And": and_, "Or": or_})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
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
        is_purine = tree.outputs.boolean("is_purine")
        is_pyrimidine = tree.outputs.boolean("is_pyrimidine")

        boolean = g.Boolean(boolean=True)
        group = ResidueName()
        index_switch = g.IndexSwitch.boolean(
            group,
            (
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                boolean,
                False,
                boolean,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                boolean,
                False,
                boolean,
            ),
        )
        BooleanAndOr(and_=and_, or_=or_, boolean=index_switch) >> is_pyrimidine
        index_switch_1 = g.IndexSwitch.boolean(
            group,
            (
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                boolean,
                False,
                boolean,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                boolean,
                False,
                boolean,
                False,
            ),
        )
        BooleanAndOr(and_=and_, or_=or_, boolean=index_switch_1) >> is_purine


ASSET = SelectNucleicType

ASSET_METADATA = {
    "catalog_id": "bd1f205b-fea5-4700-b2c2-754f3321e969",
}
