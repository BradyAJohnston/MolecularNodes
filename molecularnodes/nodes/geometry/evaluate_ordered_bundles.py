# Node-group asset "Evaluate Ordered Bundles" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BundleSocket,
    GeometrySocket,
    PackageLibrary,
    SocketAccessor,
    StringSocket,
)
from nodebpy.types import InputBundle, InputGeometry, InputString
from ._shared.sorted_bundle_paths import SortedBundlePaths


class EvaluateOrderedBundles(AssetGeometryGroup):
    """
    Evaluate Ordered Bundles

    Parameters
    ----------
    geometry : InputGeometry
        Geometry
    bundles : InputBundle
        Bundles
    prefix : InputString
        Prefix

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry
    i.bundles : BundleSocket
        Bundles
    i.prefix : StringSocket
        Prefix

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "Evaluate Ordered Bundles"
    _asset_name = "Evaluate Ordered Bundles"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""
        bundles: BundleSocket
        """Bundles"""
        prefix: StringSocket
        """Prefix"""

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
        bundles: InputBundle = None,
        prefix: InputString = "MN*",
    ):
        super().__init__(**{"Geometry": geometry, "Bundles": bundles, "Prefix": prefix})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry("Geometry")
        bundles = tree.inputs.bundle("Bundles")
        prefix = tree.inputs.string("Prefix", "MN*", optional_label=True)
        geometry_1 = tree.outputs.geometry("Geometry")

        group = SortedBundlePaths(bundle=bundles, bundle_type=prefix)
        repeat_zone = g.RepeatZone(group.o.list.list_length())
        geometry_2 = repeat_zone.items.geometry("Geometry", geometry)
        join_strings = g.JoinStrings(
            (group.o.list[repeat_zone.iteration], g.String(string="closure")),
            delimiter="/",
        )
        evaluate_closure = g.EvaluateClosure(
            g.GetBundleItem.closure(bundles, join_strings).o.item
        )
        evaluate_closure.inputs.geometry("Geometry", geometry_2.current)
        geometry_3 = evaluate_closure.outputs.geometry("Geometry")
        geometry_3 >> geometry_2.next

        geometry_2.result >> geometry_1


ASSET = EvaluateOrderedBundles

ASSET_METADATA = {
    "catalog_id": "a1e4128a-131f-4e0e-b54e-81f863aba707",
}
