# Node-group asset "Atom ID" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    AssetGeometryGroup,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputInteger
from ._shared.attribute_at_index import AttributeAtIndex


class AtomID(AssetGeometryGroup):
    """
    Atom ID

    Parameters
    ----------
    index : InputInteger
        Index

    Inputs
    ------
    i.index : IntegerSocket
        Index

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
        index: IntegerSocket
        """Index"""

    class _Outputs(SocketAccessor):
        atom_id: IntegerSocket
        """The `atom_id` attribute, read from geometry. Corresponds to the order of the atoms when they were read from the file"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        index: InputInteger = 0,
    ):
        super().__init__(**{"Index": index})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        index = tree.inputs.integer(
            "Index", 0, min_value=0, hide_value=True, default_input="INDEX"
        )
        atom_id = tree.outputs.integer(
            "atom_id",
            description="The `atom_id` attribute, read from geometry. Corresponds to the order of the atoms when they were read from the file",
        )

        AttributeAtIndex(index=index, name="atom_id") >> atom_id


ASSET = AtomID

ASSET_METADATA = {
    "catalog_id": "dfef0d3c-e718-420a-8b22-e7c3a3a9e333",
}
