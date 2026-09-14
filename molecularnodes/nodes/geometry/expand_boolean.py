# Node-group asset "Expand Boolean" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputInteger
from .offset_boolean import OffsetBoolean


class ExpandBoolean(AssetGeometryGroup):
    """
    Expand Boolean

    Parameters
    ----------
    boolean : InputBoolean
        Boolean
    expand : InputInteger
        Expand

    Inputs
    ------
    i.boolean : BooleanSocket
        Boolean
    i.expand : IntegerSocket
        Expand

    Outputs
    -------
    o.boolean : BooleanSocket
        Boolean
    """

    _name = "Expand Boolean"
    _asset_name = "Expand Boolean"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {"node_tool_idname": "geometry.expand_boolean"}

    class _Inputs(SocketAccessor):
        boolean: BooleanSocket
        """Boolean"""
        expand: IntegerSocket
        """Expand"""

    class _Outputs(SocketAccessor):
        boolean: BooleanSocket
        """Boolean"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        boolean: InputBoolean = False,
        expand: InputInteger = 0,
    ):
        super().__init__(**{"Boolean": boolean, "Expand": expand})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        boolean = tree.inputs.boolean("Boolean", False, hide_value=True)
        expand = tree.inputs.integer("Expand", 0, min_value=-2147483647)
        boolean_1 = tree.outputs.boolean("Boolean")

        boolean_math = OffsetBoolean(
            boolean=boolean, offset=expand
        ).o.boolean | OffsetBoolean(boolean=boolean, offset=-expand)
        (boolean | boolean_math) >> boolean_1


ASSET = ExpandBoolean

ASSET_METADATA = {
    "catalog_id": "0d42d20f-33f0-464a-b4c4-9fb278826421",
}
