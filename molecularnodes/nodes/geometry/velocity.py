# Node-group asset 'Velocity' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)


class Velocity(AssetGeometryGroup):
    """
    Velocity

    Outputs
    -------
    o.velocity : VectorSocket
        velocity
    """

    _name = "Velocity"
    _asset_name = "Velocity"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        velocity: VectorSocket

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self):
        super().__init__()

    def _build_group(self, tree):
        velocity = tree.outputs.vector("velocity")

        named_attribute = g.NamedAttribute.vector("velocity")

        named_attribute >> velocity


ASSET = Velocity

ASSET_METADATA = {
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
