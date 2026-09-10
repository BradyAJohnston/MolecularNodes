# Node-group asset 'VDW Radii' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputInteger
from ._shared.mn_world_scale import MN_world_scale
from .fallback_float import FallbackFloat


class VDWRadii(AssetGeometryGroup):
    """
    VDW Radii

    Parameters
    ----------
    index : InputInteger
        Index

    Inputs
    ------
    i.index : IntegerSocket
        Index

    Outputs
    -------
    o.vdw_radii : FloatSocket
        Read the `vdw_radii` attribute from the geometry
    """

    _name = "VDW Radii"
    _asset_name = "VDW Radii"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.vdw_radii"}

    class _Inputs(SocketAccessor):
        index: IntegerSocket
        """Index"""

    class _Outputs(SocketAccessor):
        vdw_radii: FloatSocket
        """Read the `vdw_radii` attribute from the geometry"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        index: InputInteger = 0,
    ):
        super().__init__(**{"Index": index})

    def _build_group(self, tree):
        index = tree.inputs.integer("Index", 0, min_value=0, default_input="INDEX")
        vdw_radii = tree.outputs.float(
            "vdw_radii", description="Read the `vdw_radii` attribute from the geometry"
        )

        (
            FallbackFloat(name="vdw_radii", fallback=MN_world_scale()).o.value.point.at(
                index
            )
            >> vdw_radii
        )


ASSET = VDWRadii

ASSET_METADATA = {
    "catalog_id": "dfef0d3c-e718-420a-8b22-e7c3a3a9e333",
}
