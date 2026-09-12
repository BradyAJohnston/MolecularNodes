# Node-group asset "MN Typed Bundles" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BundleSocket,
    ClosureSocket,
    IntegerSocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
    StringSocket,
)
from nodebpy.types import InputClosure, InputInteger, InputMenu, InputString


class MNTypedBundles(AssetGeometryGroup):
    """
    MN Typed Bundles

    Parameters
    ----------
    type : InputMenu | Literal["MN.MeshProcess"]
        Type
    closure : InputClosure
        Closure
    step : InputInteger
        step
    path : InputString
        Path

    Inputs
    ------
    i.type : MenuSocket
        Type
    i.closure : ClosureSocket
        Closure
    i.step : IntegerSocket
        step
    i.path : StringSocket
        Path

    Outputs
    -------
    o.bundle : BundleSocket
        Bundle
    """

    _name = "MN Typed Bundles"
    _asset_name = "MN Typed Bundles"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"

    class _Inputs(SocketAccessor):
        type: MenuSocket
        """Type"""
        closure: ClosureSocket
        """Closure"""
        step: IntegerSocket
        path: StringSocket
        """Path"""

    class _Outputs(SocketAccessor):
        bundle: BundleSocket
        """Bundle"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        type: InputMenu | Literal["MN.MeshProcess"] = "MN.MeshProcess",
        closure: InputClosure = None,
        step: InputInteger = 1,
        path: InputString = "",
    ):
        super().__init__(
            **{"Type": type, "Closure": closure, "step": step, "Path": path}
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        type = tree.inputs.menu("Type", optional_label=True)
        closure = tree.inputs.closure("Closure")
        step = tree.inputs.integer("step", 1)
        path = tree.inputs.string("Path", "", optional_label=True)
        bundle = tree.outputs.bundle("Bundle")

        combine_bundle = g.CombineBundle()
        combine_bundle.items.string(
            "Type",
            g.MenuSwitch.string(type, {"MN.MeshProcess": "MN.MeshProcess"}).o.output,
        )
        combine_bundle.items.closure("closure", closure)
        combine_bundle.items.integer("step", step)
        store_bundle_item = g.StoreBundleItem.bundle(
            path=path, item=combine_bundle.o.bundle
        )

        store_bundle_item >> bundle

        type.default_value = "MN.MeshProcess"


ASSET = MNTypedBundles

ASSET_METADATA = {
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
