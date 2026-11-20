# Node-group asset "Atom ID" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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


class AtomID(AssetGeometryGroup):
    """
    Atom ID

    Outputs
    -------
    o.atom_id : IntegerSocket
        The `atom_id` attribute, read from geometry. Corresponds to the order of the atoms when they were read from the file
    """

    _name = "Atom ID"
    _asset_name = "Atom ID"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.atom_id"}

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        atom_id: IntegerSocket
        """The `atom_id` attribute, read from geometry. Corresponds to the order of the atoms when they were read from the file"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self) -> None:
        super().__init__()

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atom_id = tree.outputs.integer(
            "atom_id",
            description="The `atom_id` attribute, read from geometry. Corresponds to the order of the atoms when they were read from the file",
        )

        named_attribute = g.NamedAttribute.integer("atom_id")

        named_attribute >> atom_id


ASSET = AtomID

ASSET_METADATA = {
    "catalog_id": "dfef0d3c-e718-420a-8b22-e7c3a3a9e333",
}
