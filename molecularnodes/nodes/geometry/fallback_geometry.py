# Node-group asset 'Fallback Geometry' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy.builder import (
    AssetGeometryGroup,
    GeometrySocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputGeometry
from .contains_geometry import ContainsGeometry


class FallbackGeometry(AssetGeometryGroup):
    """
    Fallback Geometry

    Parameters
    ----------
    geometry : InputGeometry
        Geometry
    fallback : InputGeometry
        Fallback

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry
    i.fallback : GeometrySocket
        Fallback

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "Fallback Geometry"
    _asset_name = "Fallback Geometry"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""
        fallback: GeometrySocket
        """Fallback"""

    class _Outputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        geometry: InputGeometry = None,
        fallback: InputGeometry = None,
    ):
        super().__init__(**{"Geometry": geometry, "Fallback": fallback})

    def _build_group(self, tree):
        geometry = tree.inputs.geometry("Geometry")
        fallback = tree.inputs.geometry("Fallback")
        geometry_1 = tree.outputs.geometry("Geometry")

        (
            ContainsGeometry(geometry=geometry).o.not_empty.switch.geometry(
                fallback, geometry
            )
            >> geometry_1
        )


ASSET = FallbackGeometry

ASSET_METADATA = {
    "catalog_id": "a1e4128a-131f-4e0e-b54e-81f863aba707",
}
