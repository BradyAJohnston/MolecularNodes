# Node-group asset 'Rotation RELION' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    PackageLibrary,
    RotationSocket,
    SocketAccessor,
)
from .tem_rotation import TEMRotation


class RotationRELION(AssetGeometryGroup):
    """
    Rotation RELION

    Outputs
    -------
    o.rotation : RotationSocket
        Rotation
    o.is_valid : BooleanSocket
        Is Valid
    """

    _name = "Rotation RELION"
    _asset_name = "Rotation RELION"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        rotation: RotationSocket
        """Rotation"""
        is_valid: BooleanSocket
        """Is Valid"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self):
        super().__init__()

    def _build_group(self, tree):
        rotation = tree.outputs.rotation("Rotation")
        is_valid = tree.outputs.boolean("Is Valid")

        group = TEMRotation()

        group >> rotation
        group.o.boolean >> is_valid


ASSET = RotationRELION

ASSET_METADATA = {
    "catalog_id": "a484cee9-1c7f-4bf8-a31c-6ffa99912ec0",
}
