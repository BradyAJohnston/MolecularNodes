# Node-group asset 'Get Geometry Atoms' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BundleSocket,
    GeometrySocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputGeometry


class GetGeometryAtoms(AssetGeometryGroup):
    """
    Get Geometry Atoms

    Parameters
    ----------
    geometry : InputGeometry
        Geometry to get the bundle of

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry to get the bundle of

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    o.bundle : BundleSocket
        Bundle
    o.atoms : GeometrySocket
        Atoms
    """

    _name = "Get Geometry Atoms"
    _asset_name = "Get Geometry Atoms"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry to get the bundle of"""

    class _Outputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""
        bundle: BundleSocket
        """Bundle"""
        atoms: GeometrySocket
        """Atoms"""

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

    def _build_group(self, tree):
        geometry = tree.inputs.geometry(
            "Geometry", description="Geometry to get the bundle of"
        )
        geometry_1 = tree.outputs.geometry("Geometry")
        bundle = tree.outputs.bundle("Bundle")
        atoms = tree.outputs.geometry("Atoms")

        get_geometry_bundle = g.GetGeometryBundle(geometry=geometry)
        get_bundle_item = g.GetBundleItem.geometry(
            get_geometry_bundle.o.bundle, g.String(string="MN/Atoms"), True
        )
        (
            get_bundle_item.o.exists.switch.geometry(
                get_geometry_bundle, get_bundle_item.o.item
            )
            >> atoms
        )
        get_bundle_item.o.exists.switch.geometry(true=get_geometry_bundle) >> geometry_1

        get_bundle_item >> bundle


ASSET = GetGeometryAtoms

ASSET_METADATA = {
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
