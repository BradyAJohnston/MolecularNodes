# Node-group asset "oxDNA Offset" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
    VectorSocket,
)
from nodebpy.types import InputInteger
from .angstrom_to_world import AngstromToWorld
from .oxdna_normal import OxDNANormal
from .oxdna_vector import OxDNAVector


class OxDNAOffset(AssetGeometryGroup):
    """
    oxDNA Offset

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
    o.offset : VectorSocket
        Offset
    """

    _name = "oxDNA Offset"
    _asset_name = "oxDNA Offset"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        index: IntegerSocket
        """Index"""

    class _Outputs(SocketAccessor):
        offset: VectorSocket
        """Offset"""

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
        offset = tree.outputs.vector("Offset")

        group = OxDNAVector(index=index)
        (
            group.o.base_vector * -0.34
            + group.o.base_vector.cross(OxDNANormal(index=index))
            * AngstromToWorld(angstrom=3.408)
            >> offset
        )


ASSET = OxDNAOffset

ASSET_METADATA = {
    "catalog_id": "0094c3e0-7885-427b-81b4-187a84dcff18",
}
