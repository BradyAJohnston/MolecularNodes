# Node-group asset "Select Attribute" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    PackageLibrary,
    SocketAccessor,
    StringSocket,
)
from nodebpy.types import InputBoolean, InputString
from ._shared.mn_select_attribute import MN_select_attribute


class SelectAttribute(AssetGeometryGroup):
    """
    Select Attribute

    Parameters
    ----------
    and_ : InputBoolean
        The resulting selection must overlap with this input selection
    or_ : InputBoolean
        The resulting selection can be calculated from this node or be from this input selection
    name : InputString
        Name

    Inputs
    ------
    i.and_ : BooleanSocket
        The resulting selection must overlap with this input selection
    i.or_ : BooleanSocket
        The resulting selection can be calculated from this node or be from this input selection
    i.name : StringSocket
        Name

    Outputs
    -------
    o.selection : BooleanSocket
        The calculated selection
    o.inverted : BooleanSocket
        The inverse of the calculated selection
    """

    _name = "Select Attribute"
    _asset_name = "Select Attribute"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.select_attribute"}

    class _Inputs(SocketAccessor):
        and_: BooleanSocket
        """The resulting selection must overlap with this input selection"""
        or_: BooleanSocket
        """The resulting selection can be calculated from this node or be from this input selection"""
        name: StringSocket
        """Name"""

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
        name: InputString = "is_peptide",
    ):
        super().__init__(**{"And": and_, "Or": or_, "Name": name})

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
        name = tree.inputs.string("Name", "is_peptide", optional_label=True)
        selection = tree.outputs.boolean(
            "Selection", description="The calculated selection"
        )
        inverted = tree.outputs.boolean(
            "Inverted", description="The inverse of the calculated selection"
        )

        group = MN_select_attribute(and_=and_, or_=or_, name=name)

        group >> selection
        group.o.inverted >> inverted


ASSET = SelectAttribute

ASSET_METADATA = {
    "catalog_id": "bd1f205b-fea5-4700-b2c2-754f3321e969",
}
