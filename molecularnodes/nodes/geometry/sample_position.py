# Node-group asset 'Sample Position' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    GeometrySocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputGeometry, InputInteger, InputVector


class SamplePosition(AssetGeometryGroup):
    """
    A convenience wrapper around the `Sample Index` and `Position` nodes

    Parameters
    ----------
    geometry : InputGeometry
        The geometry to sample the `Position` from
    position : InputVector
        The `Position` field to sample the values from
    index : InputInteger
        The `Index` at which to sample the `Position` field from

    Inputs
    ------
    i.geometry : GeometrySocket
        The geometry to sample the `Position` from
    i.position : VectorSocket
        The `Position` field to sample the values from
    i.index : IntegerSocket
        The `Index` at which to sample the `Position` field from

    Outputs
    -------
    o.position : VectorSocket
        The sampled `Position` field
    """

    _name = "Sample Position"
    _asset_name = "Sample Position"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "description": "A convenience wrapper around the `Sample Index` and `Position` nodes",
        "node_tool_idname": "geometry.sample_position",
    }

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """The geometry to sample the `Position` from"""
        position: VectorSocket
        """The `Position` field to sample the values from"""
        index: IntegerSocket
        """The `Index` at which to sample the `Position` field from"""

    class _Outputs(SocketAccessor):
        position: VectorSocket
        """The sampled `Position` field"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        geometry: InputGeometry = None,
        position: InputVector = None,
        index: InputInteger = 0,
    ):
        super().__init__(**{"Geometry": geometry, "Position": position, "Index": index})

    def _build_group(self, tree):
        geometry = tree.inputs.geometry(
            "Geometry", description="The geometry to sample the `Position` from"
        )
        position = tree.inputs.vector(
            "Position",
            (0.0, 0.0, 0.0),
            description="The `Position` field to sample the values from",
            hide_value=True,
            default_input="POSITION",
        )
        index = tree.inputs.integer(
            "Index",
            0,
            description="The `Index` at which to sample the `Position` field from",
            default_input="INDEX",
        )
        position_1 = tree.outputs.vector(
            "Position", description="The sampled `Position` field"
        )

        (
            geometry
            >> g.SampleIndex(
                value=position, index=index, data_type="FLOAT_VECTOR", clamp=True
            )
            >> position_1
        )


ASSET = SamplePosition

ASSET_METADATA = {
    "description": "A convenience wrapper around the `Sample Index` and `Position` nodes",
    "catalog_id": "dd5f0199-fa8b-4b01-a972-2dc586a3e60f",
}
