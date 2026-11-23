# Node-group asset "Atom Name" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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


class AtomName(AssetGeometryGroup):
    """
    Geometry Nodes doesn't currently support text attributes, so strings like atom names have to be first converted to integers and mapping to atom names back and forth. Definitions for the atom names are available on the GitHub page

    Outputs
    -------
    o.atom_name : IntegerSocket
        The `atom_name` attribute read from the points, an integer representation of the atom names. Corresponds to `CA` and `O` in the file etc.
    """

    _name = "Atom Name"
    _asset_name = "Atom Name"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"
    _tree_properties = {
        "description": "Geometry Nodes doesn't currently support text attributes, so strings like atom names have to be first converted to integers and mapping to atom names back and forth. Definitions for the atom names are available on the GitHub page",
        "node_tool_idname": "geometry.atom_name",
    }

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        atom_name: IntegerSocket
        """The `atom_name` attribute read from the points, an integer representation of the atom names. Corresponds to `CA` and `O` in the file etc."""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self) -> None:
        super().__init__()

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        atom_name = tree.outputs.integer(
            "atom_name",
            description="The `atom_name` attribute read from the points, an integer representation of the atom names. Corresponds to `CA` and `O` in the file etc. ",
        )

        named_attribute = g.NamedAttribute.integer("atom_name")

        named_attribute >> atom_name


ASSET = AtomName

ASSET_METADATA = {
    "description": "Geometry Nodes doesn't currently support text attributes, so strings like atom names have to be first converted to integers and mapping to atom names back and forth. Definitions for the atom names are available on the GitHub page",
    "catalog_id": "dfef0d3c-e718-420a-8b22-e7c3a3a9e333",
}
