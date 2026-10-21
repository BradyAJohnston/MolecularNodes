# Node-group asset "Set Geometry Volume" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
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


class SetGeometryVolume(AssetGeometryGroup):
    """
    Set Geometry Volume

    Parameters
    ----------
    geometry : InputGeometry
        Styled geometry, such as a density's iso surface
    volume : InputGeometry
        The source volume the geometry was made from

    Inputs
    ------
    i.geometry : GeometrySocket
        Styled geometry, such as a density's iso surface
    i.volume : GeometrySocket
        The source volume the geometry was made from

    Outputs
    -------
    o.geometry : GeometrySocket
        The geometry, with the source volume in its `MN/Volume` bundle
    """

    _name = "Set Geometry Volume"
    _asset_name = "Set Geometry Volume"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Styled geometry, such as a density's iso surface"""
        volume: GeometrySocket
        """The source volume the geometry was made from"""

    class _Outputs(SocketAccessor):
        geometry: GeometrySocket
        """The geometry, with the source volume in its `MN/Volume` bundle"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        geometry: InputGeometry = None,
        volume: InputGeometry = None,
    ):
        super().__init__(**{"Geometry": geometry, "Volume": volume})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry(
            "Geometry", description="Styled geometry, such as a density's iso surface"
        )
        volume = tree.inputs.geometry(
            "Volume", description="The source volume the geometry was made from"
        )
        geometry_1 = tree.outputs.geometry(
            "Geometry",
            description="The geometry, with the source volume in its `MN/Volume` bundle",
        )

        with g.Frame("Keep the source volume in the 'MN/Volume' bundle"):
            get_geometry_bundle = geometry.bundle()
            get_bundle_item = g.GetBundleItem.geometry(get_geometry_bundle, "MN/Volume")
            get_bundle_item_1 = g.GetBundleItem.geometry(volume.bundle(), "MN/Volume")
            _string = g.String(
                "Density styles turn the volume into a mesh. The volume is kept in the 'MN/Volume' bundle of the output, which is not rendered, so nodes such as Color Sample Grid can still sample it through Object Info. A volume already in the geometry's or the volume's 'MN/Volume' bundle is kept instead of being replaced."
            )
            switch = get_bundle_item.o.exists.switch.geometry(
                get_bundle_item_1.o.exists.switch.geometry(
                    volume, get_bundle_item_1.o.item
                ),
                get_bundle_item.o.item,
            )
            store_bundle_item = g.StoreBundleItem.geometry(
                get_geometry_bundle, "MN/Volume", switch
            )
        geometry.set_bundle(store_bundle_item) >> geometry_1


ASSET = SetGeometryVolume

ASSET_METADATA = {
    "description": "Keep the source volume of a density in the `MN/Volume` bundle of its styled geometry, so it can still be sampled",
    "catalog_id": "35b7cc56-45fd-4113-8e88-3c7387edb2d8",
}
