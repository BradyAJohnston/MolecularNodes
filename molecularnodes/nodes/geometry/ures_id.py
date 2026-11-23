# Node-group asset "URes ID" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)


class UResID(AssetGeometryGroup):
    """
    URes ID

    Outputs
    -------
    o.ures_id : IntegerSocket
        ures_id
    o.size : IntegerSocket
        The total of all of the values in the corresponding group
    """

    _name = "URes ID"
    _asset_name = "URes ID"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        ures_id: IntegerSocket
        size: IntegerSocket
        """The total of all of the values in the corresponding group"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self) -> None:
        super().__init__()

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        ures_id = tree.outputs.integer("ures_id")
        size = tree.outputs.integer(
            "Size",
            description="The total of all of the values in the corresponding group",
        )

        named_attribute = g.NamedAttribute.integer("ures_id")
        accumulate_field = g.AccumulateField.point.integer(
            group_id=named_attribute.o.attribute
        )

        named_attribute >> ures_id
        accumulate_field.o.total >> size


ASSET = UResID

ASSET_METADATA = {
    "catalog_id": "dfef0d3c-e718-420a-8b22-e7c3a3a9e333",
}
