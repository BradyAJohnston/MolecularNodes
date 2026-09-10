# Node-group asset 'Edge Length' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    PackageLibrary,
    SocketAccessor,
)


class EdgeLength(AssetGeometryGroup):
    """
    Edge Length

    Outputs
    -------
    o.length : FloatSocket
        Length
    """

    _name = "Edge Length"
    _asset_name = "Edge Length"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        length: FloatSocket
        """Length"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self):
        super().__init__()

    def _build_group(self, tree):
        length = tree.outputs.float("Length")

        edge_vertices = g.EdgeVertices()
        (
            edge_vertices.o.position_1.distance(
                edge_vertices.o.position_2
            ).edge.evaluate()
            >> length
        )


ASSET = EdgeLength

ASSET_METADATA = {
    "description": "Returns the length of an edge",
    "copyright": "Blender Foundation",
    "license": "CC0 - Public Domain",
    "catalog_id": "b6bb38bb-bfe1-4a12-b2bc-fc87d060864f",
}
