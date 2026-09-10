# Node-group asset 'oxDNA Rotation' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    IntegerSocket,
    PackageLibrary,
    RotationSocket,
    SocketAccessor,
)
from nodebpy.types import InputInteger
from .oxdna_normal import OxDNANormal
from .oxdna_vector import OxDNAVector


class OxDNARotation(AssetGeometryGroup):
    """
    oxDNA Rotation

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
    o.rotation : RotationSocket
        Rotation
    """

    _name = "oxDNA Rotation"
    _asset_name = "oxDNA Rotation"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        index: IntegerSocket
        """Index"""

    class _Outputs(SocketAccessor):
        rotation: RotationSocket
        """Rotation"""

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
        index = tree.inputs.integer("Index", 0, min_value=0, default_input="INDEX")
        rotation = tree.outputs.rotation("Rotation")

        axes_to_rotation = g.AxesToRotation(
            primary_axis=OxDNAVector(index=index),
            secondary_axis=OxDNANormal(index=index),
        )

        axes_to_rotation >> rotation


ASSET = OxDNARotation

ASSET_METADATA = {
    "catalog_id": "0094c3e0-7885-427b-81b4-187a84dcff18",
}
