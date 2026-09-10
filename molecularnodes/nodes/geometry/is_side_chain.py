# Node-group asset "Is Side Chain" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from ._shared.mn_select_nucleic import MN_select_nucleic
from ._shared.mn_select_peptide import MN_select_peptide
from .boolean_andor import BooleanAndOr
from .fallback_boolean import FallbackBoolean
from .is_alpha_carbon import IsAlphaCarbon


class IsSideChain(AssetGeometryGroup):
    """
    Is Side Chain

    Parameters
    ----------
    and_ : InputBoolean
        The resulting selection must overlap with this input selection
    or_ : InputBoolean
        The resulting selection can be calculated from this node or be from this input selection
    include_ca : InputBoolean
        Include the alpha carbon as part of the side chain

    Inputs
    ------
    i.and_ : BooleanSocket
        The resulting selection must overlap with this input selection
    i.or_ : BooleanSocket
        The resulting selection can be calculated from this node or be from this input selection
    i.include_ca : BooleanSocket
        Include the alpha carbon as part of the side chain

    Outputs
    -------
    o.selection : BooleanSocket
        The calculated selection
    o.inverted : BooleanSocket
        The inverse of the calculated selection
    """

    _name = "Is Side Chain"
    _asset_name = "Is Side Chain"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.is_side_chain"}

    class _Inputs(SocketAccessor):
        and_: BooleanSocket
        """The resulting selection must overlap with this input selection"""
        or_: BooleanSocket
        """The resulting selection can be calculated from this node or be from this input selection"""
        include_ca: BooleanSocket
        """Include the alpha carbon as part of the side chain"""

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
        include_ca: InputBoolean = True,
    ):
        super().__init__(**{"And": and_, "Or": or_, "Include CA": include_ca})

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
        include_ca = tree.inputs.boolean(
            "Include CA",
            True,
            description="Include the alpha carbon as part of the side chain",
        )
        selection = tree.outputs.boolean(
            "Selection", description="The calculated selection"
        )
        inverted = tree.outputs.boolean(
            "Inverted", description="The inverse of the calculated selection"
        )

        group = FallbackBoolean(
            name="is_side_chain",
            fallback=MN_select_nucleic().o.is_side_chain
            | MN_select_peptide().o.is_side_chain,
        )
        switch = include_ca.switch.boolean(
            g.BooleanMath.subtract(group, IsAlphaCarbon().o.selection), group
        )
        group_1 = BooleanAndOr(and_=and_, or_=or_, boolean=switch)

        group_1 >> selection
        group_1.o.inverted >> inverted


ASSET = IsSideChain

ASSET_METADATA = {
    "catalog_id": "bd1f205b-fea5-4700-b2c2-754f3321e969",
}
