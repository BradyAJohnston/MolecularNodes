# Node-group asset 'Sample Mix Vector' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    GeometrySocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputFloat, InputGeometry, InputInteger, InputVector


class SampleMixVector(AssetGeometryGroup):
    """
    Sample Mix Vector

    Parameters
    ----------
    a : InputGeometry
        Geometry A to sample and mix from
    b : InputGeometry
        Geometry B to sample and mix to
    position : InputVector
        The field to sample from each geometry, defaulting to `Position`
    factor : InputFloat
        The amount to mix from A to B
    index : InputInteger
        `Index` on the geometries to sample from

    Inputs
    ------
    i.a : GeometrySocket
        Geometry A to sample and mix from
    i.b : GeometrySocket
        Geometry B to sample and mix to
    i.position : VectorSocket
        The field to sample from each geometry, defaulting to `Position`
    i.factor : FloatSocket
        The amount to mix from A to B
    i.index : IntegerSocket
        `Index` on the geometries to sample from

    Outputs
    -------
    o.vector : VectorSocket
        The final mixed vector
    """

    _name = "Sample Mix Vector"
    _asset_name = "Sample Mix Vector"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry.sample_mix_vector"}

    class _Inputs(SocketAccessor):
        a: GeometrySocket
        """Geometry A to sample and mix from"""
        b: GeometrySocket
        """Geometry B to sample and mix to"""
        position: VectorSocket
        """The field to sample from each geometry, defaulting to `Position`"""
        factor: FloatSocket
        """The amount to mix from A to B"""
        index: IntegerSocket
        """`Index` on the geometries to sample from"""

    class _Outputs(SocketAccessor):
        vector: VectorSocket
        """The final mixed vector"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        a: InputGeometry = None,
        b: InputGeometry = None,
        position: InputVector = None,
        factor: InputFloat = 0.5,
        index: InputInteger = 0,
    ):
        super().__init__(
            **{"A": a, "B": b, "Position": position, "Factor": factor, "Index": index}
        )

    def _build_group(self, tree):
        a = tree.inputs.geometry("A", description="Geometry A to sample and mix from")
        b = tree.inputs.geometry("B", description="Geometry B to sample and mix to")
        position = tree.inputs.vector(
            "Position",
            (0.0, 0.0, 0.0),
            description="The field to sample from each geometry, defaulting to `Position`",
            hide_value=True,
            default_input="POSITION",
        )
        factor = tree.inputs.float(
            "Factor",
            0.5,
            description="The amount to mix from A to B",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        index = tree.inputs.integer(
            "Index",
            0,
            description="`Index` on the geometries to sample from",
            default_input="INDEX",
        )
        vector = tree.outputs.vector("Vector", description="The final mixed vector")

        mix = g.Mix(
            factor_float=factor,
            a_vector=g.SampleIndex(
                geometry=a, value=position, index=index, data_type="FLOAT_VECTOR"
            ),
            b_vector=g.SampleIndex(
                geometry=b, value=position, index=index, data_type="FLOAT_VECTOR"
            ),
            data_type="VECTOR",
            clamp_factor=True,
        )

        mix.o.result_vector >> vector


ASSET = SampleMixVector

ASSET_METADATA = {
    "catalog_id": "dd5f0199-fa8b-4b01-a972-2dc586a3e60f",
}
