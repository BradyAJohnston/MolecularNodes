# Node-group asset "Check Geometry" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    GeometrySocket,
    PackageLibrary,
    SocketAccessor,
    StringSocket,
)
from nodebpy.types import InputGeometry, InputString
from .contains_geometry import ContainsGeometry


class CheckGeometry(AssetGeometryGroup):
    """
    Check Geometry

    Parameters
    ----------
    geometry : InputGeometry
        Geometry
    message : InputString
        Message

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry
    i.message : StringSocket
        Message

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "Check Geometry"
    _asset_name = "Check Geometry"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""
        message: StringSocket
        """Message"""

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
        message: InputString = "Input contains no geometry, check your node connections.",
    ):
        super().__init__(**{"Geometry": geometry, "Message": message})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry("Geometry")
        message = tree.inputs.string(
            "Message",
            "Input contains no geometry, check your node connections.",
            optional_label=True,
        )
        geometry_1 = tree.outputs.geometry("Geometry")

        _warning = g.Warning(
            show=ContainsGeometry(geometry=geometry).o.empty, message=message
        )

        geometry >> geometry_1


ASSET = CheckGeometry

ASSET_METADATA = {
    "catalog_id": "a1e4128a-131f-4e0e-b54e-81f863aba707",
}
