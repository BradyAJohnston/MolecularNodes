# Node-group asset 'Geoemtry to Planar' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    GeometrySocket,
    MatrixSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputGeometry
from ._shared.geometry_principal_components import GeometryPrincipalComponents


class GeoemtryToPlanar(AssetGeometryGroup):
    """
    Geoemtry to Planar

    Parameters
    ----------
    geometry : InputGeometry
        Geometry to transform
    selection : InputBoolean
        The parts of the geometry that contibute to the planar calculation

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry to transform
    i.selection : BooleanSocket
        The parts of the geometry that contibute to the planar calculation

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    o.transform : MatrixSocket
        Transform
    """

    _name = "Geoemtry to Planar"
    _asset_name = "Geoemtry to Planar"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry to transform"""
        selection: BooleanSocket
        """The parts of the geometry that contibute to the planar calculation"""

    class _Outputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""
        transform: MatrixSocket
        """Transform"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        geometry: InputGeometry = None,
        selection: InputBoolean = True,
    ):
        super().__init__(**{"Geometry": geometry, "Selection": selection})

    def _build_group(self, tree):
        geometry = tree.inputs.geometry("Geometry", description="Geometry to transform")
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="The parts of the geometry that contibute to the planar calculation",
            hide_value=True,
        )
        geometry_1 = tree.outputs.geometry("Geometry")
        transform = tree.outputs.matrix("Transform")

        group = GeometryPrincipalComponents(
            geometry=g.SeparateGeometry.point(geometry, selection).o.selection
        )
        with g.Frame("Transform to Planar"):
            combine_transform = g.CombineTransform(
                translation=group.o.group_center * -1.0,
                rotation=group.o.rotation.invert(),
            )
        transform_geometry = g.TransformGeometry(
            geometry=geometry, transform=combine_transform, mode="Matrix"
        )

        transform_geometry >> geometry_1
        combine_transform >> transform


ASSET = GeoemtryToPlanar

ASSET_METADATA = {
    "catalog_id": "a1e4128a-131f-4e0e-b54e-81f863aba707",
}
