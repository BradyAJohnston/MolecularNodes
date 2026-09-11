# Node-group asset "Clear Instance Transforms" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    GeometrySocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputGeometry


class ClearInstanceTransforms(AssetGeometryGroup):
    """
    Clear Instance Transforms

    Parameters
    ----------
    geometry : InputGeometry
        Geometry
    selection : InputBoolean
        Selection
    realize_all : InputBoolean
        Realize all levels of nested instances for a top-level instances. Overrides the value of the Depth input

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry
    i.selection : BooleanSocket
        Selection
    i.realize_all : BooleanSocket
        Realize all levels of nested instances for a top-level instances. Overrides the value of the Depth input

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "Clear Instance Transforms"
    _asset_name = "Clear Instance Transforms"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""
        selection: BooleanSocket
        """Selection"""
        realize_all: BooleanSocket
        """Realize all levels of nested instances for a top-level instances. Overrides the value of the Depth input"""

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
        selection: InputBoolean = True,
        realize_all: InputBoolean = False,
    ):
        super().__init__(
            **{"Geometry": geometry, "Selection": selection, "Realize All": realize_all}
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry("Geometry")
        selection = tree.inputs.boolean("Selection", True, hide_value=True)
        realize_all = tree.inputs.boolean(
            "Realize All",
            False,
            description="Realize all levels of nested instances for a top-level instances. Overrides the value of the Depth input",
        )
        geometry_1 = tree.outputs.geometry("Geometry")

        (
            geometry
            >> g.SetInstanceTransform(selection=selection, transform=g.CombineMatrix())
            >> g.RealizeInstances(realize_all=realize_all)
            >> geometry_1
        )


ASSET = ClearInstanceTransforms

ASSET_METADATA = {
    "catalog_id": "a1e4128a-131f-4e0e-b54e-81f863aba707",
}
