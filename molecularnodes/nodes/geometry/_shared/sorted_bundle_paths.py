# Node group "Sorted Bundle Paths" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    BundleSocket,
    CustomGeometryGroup,
    SocketAccessor,
    StringSocket,
)
from nodebpy.types import InputBundle, InputString


class SortedBundlePaths(CustomGeometryGroup):
    """
    Sorted Bundle Paths

    Parameters
    ----------
    bundle : InputBundle
        Bundle
    bundle_type : InputString
        Bundle Type

    Inputs
    ------
    i.bundle : BundleSocket
        Bundle
    i.bundle_type : StringSocket
        Bundle Type

    Outputs
    -------
    o.list : StringSocket
        List
    """

    _name = "Sorted Bundle Paths"

    class _Inputs(SocketAccessor):
        bundle: BundleSocket
        """Bundle"""
        bundle_type: StringSocket
        """Bundle Type"""

    class _Outputs(SocketAccessor):
        list: StringSocket
        """List"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        bundle: InputBundle = None,
        bundle_type: InputString = "MN*",
    ):
        super().__init__(**{"Bundle": bundle, "Bundle Type": bundle_type})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        bundle = tree.inputs.bundle("Bundle")
        bundle_type = tree.inputs.string("Bundle Type", "MN*", optional_label=True)
        list = tree.outputs.string("List")

        closure_zone = g.ClosureZone()
        index = closure_zone.inputs.integer("Index", structure_type="SINGLE")
        item = closure_zone.outputs.integer("Item", structure_type="DYNAMIC")
        get_nested_bundle_paths = g.GetNestedBundlePaths(
            bundle=bundle,
            bundle_type=bundle_type,
            mode="Bundle Type",
            pattern_mode="Wildcard",
        )
        format_string = g.FormatString("{v}/step", items={"v": get_nested_bundle_paths})
        g.GetBundleItem.integer(bundle, format_string.o.string[index]).o.item >> item
        closure_to_list = g.ClosureToList(
            count=format_string.o.string.list_length(), closure=closure_zone.closure
        )
        item_1 = closure_to_list.items.integer("Item", structure_type="DYNAMIC")
        get_nested_bundle_paths.o.paths.sort(item_1) >> list
