# Node-group asset 'Atom Name' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy.builder import (
    AssetGeometryGroup,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputInteger
from ._shared.attribute_at_index import AttributeAtIndex


class AtomName(AssetGeometryGroup):
    """
    Geometry Nodes doesn't currently support text attributes, so strings like atom names have to be first converted to integers and mapping to atom names back and forth. Definitions for the atom names are available on the GitHub page

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
    o.atom_name : IntegerSocket
        The `atom_name` attribute read from the points, an integer representation of the atom names. Corresponds to `CA` and `O` in the file etc.
    """

    _name = "Atom Name"
    _asset_name = "Atom Name"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {
        "description": "Geometry Nodes doesn't currently support text attributes, so strings like atom names have to be first converted to integers and mapping to atom names back and forth. Definitions for the atom names are available on the GitHub page",
        "node_tool_idname": "geometry.atom_name",
    }

    class _Inputs(SocketAccessor):
        index: IntegerSocket
        """Index"""

    class _Outputs(SocketAccessor):
        atom_name: IntegerSocket
        """The `atom_name` attribute read from the points, an integer representation of the atom names. Corresponds to `CA` and `O` in the file etc."""

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

    def _build_group(self, tree):
        index = tree.inputs.integer(
            "Index", 0, min_value=0, hide_value=True, default_input="INDEX"
        )
        atom_name = tree.outputs.integer(
            "atom_name",
            description="The `atom_name` attribute read from the points, an integer representation of the atom names. Corresponds to `CA` and `O` in the file etc. ",
        )

        AttributeAtIndex(index=index, name="atom_name") >> atom_name


ASSET = AtomName

ASSET_METADATA = {
    "description": "Geometry Nodes doesn't currently support text attributes, so strings like atom names have to be first converted to integers and mapping to atom names back and forth. Definitions for the atom names are available on the GitHub page",
    "catalog_id": "dfef0d3c-e718-420a-8b22-e7c3a3a9e333",
}
