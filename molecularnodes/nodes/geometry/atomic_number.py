# Node-group asset "Atomic Number" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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


class AtomicNumber(AssetGeometryGroup):
    """
    Atomic Number

    Outputs
    -------
    o.atomic_number : IntegerSocket
        Read the `atomic_number` attribute from the geometry corresponding to elements
    """

    _name = "Atomic Number"
    _asset_name = "Atomic Number"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.atomic_number"}

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        atomic_number: IntegerSocket
        """Read the `atomic_number` attribute from the geometry corresponding to elements"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self) -> None:
        super().__init__()

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atomic_number = tree.outputs.integer(
            "atomic_number",
            description="Read the `atomic_number` attribute from the geometry corresponding to elements",
        )

        named_attribute = g.NamedAttribute.integer("atomic_number")

        named_attribute >> atomic_number


ASSET = AtomicNumber

ASSET_METADATA = {
    "catalog_id": "dfef0d3c-e718-420a-8b22-e7c3a3a9e333",
}
