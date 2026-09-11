# Node-group asset "Set URes ID" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
)
from nodebpy.types import InputGeometry
from .unique_residue_id import UniqueResidueID


class SetUResID(AssetGeometryGroup):
    """
    Set URes ID

    Parameters
    ----------
    geometry : InputGeometry
        Geometry

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "Set URes ID"
    _asset_name = "Set URes ID"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""

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
    ):
        super().__init__(**{"Geometry": geometry})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry("Geometry")
        geometry_1 = tree.outputs.geometry("Geometry")

        (
            geometry
            >> g.StoreNamedAttribute.point.integer(
                name="ures_id", value=UniqueResidueID()
            )
            >> geometry_1
        )


ASSET = SetUResID

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
