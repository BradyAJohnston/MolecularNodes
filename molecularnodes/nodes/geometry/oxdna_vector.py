# Node-group asset "oxDNA Vector" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
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
    VectorSocket,
)
from nodebpy.types import InputInteger


class OxDNAVector(AssetGeometryGroup):
    """
    oxDNA Vector

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
    o.base_vector : VectorSocket
        base_vector
    """

    _name = "oxDNA Vector"
    _asset_name = "oxDNA Vector"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        index: IntegerSocket
        """Index"""

    class _Outputs(SocketAccessor):
        base_vector: VectorSocket

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
        index = tree.inputs.integer("Index", 0, min_value=0, default_input="INDEX")
        base_vector = tree.outputs.vector("base_vector")

        (
            g.NamedAttribute.vector("base_vector").o.attribute.point.at(index)
            >> base_vector
        )


ASSET = OxDNAVector

ASSET_METADATA = {
    "catalog_id": "0094c3e0-7885-427b-81b4-187a84dcff18",
}
