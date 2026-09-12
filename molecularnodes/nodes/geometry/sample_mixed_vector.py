# Node-group asset "Sample Mixed Vector" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    GeometrySocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputFloat, InputGeometry, InputVector
from .index_mix_vector import IndexMixVector


class SampleMixedVector(AssetGeometryGroup):
    """
    Sample Mixed Vector

    Parameters
    ----------
    geometry : InputGeometry
        The geometry to sample the values from
    vector : InputVector
        The field to mix and evaluate on the sample geometry
    index : InputFloat
        The index to sample the value from. The fractional component of the index is used to mix between values using `Index Mix ...` nodes

    Inputs
    ------
    i.geometry : GeometrySocket
        The geometry to sample the values from
    i.vector : VectorSocket
        The field to mix and evaluate on the sample geometry
    i.index : FloatSocket
        The index to sample the value from. The fractional component of the index is used to mix between values using `Index Mix ...` nodes

    Outputs
    -------
    o.vector : VectorSocket
        The evaluated and mixed field, sampled from the sample geometry at the given `Index`
    """

    _name = "Sample Mixed Vector"
    _asset_name = "Sample Mixed Vector"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry.sample_mixed_vector"}

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """The geometry to sample the values from"""
        vector: VectorSocket
        """The field to mix and evaluate on the sample geometry"""
        index: FloatSocket
        """The index to sample the value from. The fractional component of the index is used to mix between values using `Index Mix ...` nodes"""

    class _Outputs(SocketAccessor):
        vector: VectorSocket
        """The evaluated and mixed field, sampled from the sample geometry at the given `Index`"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        geometry: InputGeometry = None,
        vector: InputVector = None,
        index: InputFloat = 0.0,
    ):
        super().__init__(**{"Geometry": geometry, "Vector": vector, "Index": index})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry(
            "Geometry", description="The geometry to sample the values from"
        )
        vector = tree.inputs.vector(
            "Vector",
            (0.0, 0.0, 0.0),
            description="The field to mix and evaluate on the sample geometry",
            hide_value=True,
            default_input="POSITION",
        )
        index = tree.inputs.float(
            "Index",
            0.0,
            description="The index to sample the value from. The fractional component of the index is used to mix between values using `Index Mix ...` nodes",
        )
        vector_1 = tree.outputs.vector(
            "Vector",
            description="The evaluated and mixed field, sampled from the sample geometry at the given `Index`",
        )

        group = IndexMixVector(value=vector, index=index)
        (
            geometry
            >> g.SampleIndex(
                value=group.o.value,
                index=group.o.from_,
                data_type="FLOAT_VECTOR",
                clamp=True,
            )
            >> vector_1
        )


ASSET = SampleMixedVector

ASSET_METADATA = {
    "catalog_id": "dd5f0199-fa8b-4b01-a972-2dc586a3e60f",
}
