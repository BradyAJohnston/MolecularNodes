# Node-group asset "Cumulative Length" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    PackageLibrary,
    SocketAccessor,
)


class CumulativeLength(AssetGeometryGroup):
    """
    Cumulative Length

    Outputs
    -------
    o.length : FloatSocket
        The length along the current spline added to all previous spline lengths
    """

    _name = "Cumulative Length"
    _asset_name = "Cumulative Length"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.cumulative_length"}

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        length: FloatSocket
        """The length along the current spline added to all previous spline lengths"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self):
        super().__init__()

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        length = tree.outputs.float(
            "Length",
            description="The length along the current spline added to all previous spline lengths",
        )

        spline_parameter = g.SplineParameter()
        (
            (
                spline_parameter.o.length.spline.trailing() + spline_parameter.o.length
            ).point.evaluate()
            >> length
        )


ASSET = CumulativeLength

ASSET_METADATA = {
    "catalog_id": "9c167a5c-d0a6-457d-9e9c-90f3edd29e10",
}
