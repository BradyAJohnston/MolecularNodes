# Node-group asset "Attribute Run" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
    StringSocket,
)
from nodebpy.types import InputInteger, InputString
from .integer_run import IntegerRun


class AttributeRun(AssetGeometryGroup):
    """
    Group mask increments whenever the attribute or the Group ID changes

    Parameters
    ----------
    name : InputString
        The `Named Attribute` to read from the geometry before applying the `Integer Run` node
    group_id : InputInteger
        The `Group ID` to use for the `Integer Run` node

    Inputs
    ------
    i.name : StringSocket
        The `Named Attribute` to read from the geometry before applying the `Integer Run` node
    i.group_id : IntegerSocket
        The `Group ID` to use for the `Integer Run` node

    Outputs
    -------
    o.is_different : BooleanSocket
        The current point's `Named Attribute` value is different from the previous
    o.group_id : IntegerSocket
        The new `Group ID`, increasing whenever the attribute or the `Group ID` values change
    """

    _name = "Attribute Run"
    _asset_name = "Attribute Run"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {
        "description": "Group mask increments whenever the attribute or the Group ID changes",
        "node_tool_idname": "geometry.attribute_run",
    }

    class _Inputs(SocketAccessor):
        name: StringSocket
        """The `Named Attribute` to read from the geometry before applying the `Integer Run` node"""
        group_id: IntegerSocket
        """The `Group ID` to use for the `Integer Run` node"""

    class _Outputs(SocketAccessor):
        is_different: BooleanSocket
        """The current point's `Named Attribute` value is different from the previous"""
        group_id: IntegerSocket
        """The new `Group ID`, increasing whenever the attribute or the `Group ID` values change"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        name: InputString = "",
        group_id: InputInteger = 0,
    ):
        super().__init__(**{"Name": name, "Group ID": group_id})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        name = tree.inputs.string(
            "Name",
            "",
            description="The `Named Attribute` to read from the geometry before applying the `Integer Run` node",
            optional_label=True,
        )
        group_id = tree.inputs.integer(
            "Group ID",
            0,
            description="The `Group ID` to use for the `Integer Run` node",
            hide_value=True,
        )
        is_different = tree.outputs.boolean(
            "Is Different",
            description="The current point's `Named Attribute` value is different from the previous",
        )
        group_id_1 = tree.outputs.integer(
            "Group ID",
            description="The new `Group ID`, increasing whenever the attribute or the `Group ID` values change",
        )

        group = IntegerRun(
            value=g.NamedAttribute.integer(name).o.attribute, group_id=group_id
        )

        group >> is_different
        group.o.group_id >> group_id_1


ASSET = AttributeRun

ASSET_METADATA = {
    "description": "Group mask increments whenever the attribute or the Group ID changes",
    "catalog_id": "7ccb8802-a69f-483e-bf6e-4a47aaa9e940",
}
