# Node-group asset "Sample Mix Float" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputFloat, InputGeometry, InputInteger


class SampleMixFloat(AssetGeometryGroup):
    """
    Sample a float value from two different geometries and mix from A to B

    Parameters
    ----------
    a : InputGeometry
        Geometry A to sample and mix from
    b : InputGeometry
        Geometry B to sample and mix to
    value : InputFloat
        Float field to sample and mix
    factor : InputFloat
        Amount to mix from A to B
    index : InputInteger
        `Index` on the geometries to sample from

    Inputs
    ------
    i.a : GeometrySocket
        Geometry A to sample and mix from
    i.b : GeometrySocket
        Geometry B to sample and mix to
    i.value : FloatSocket
        Float field to sample and mix
    i.factor : FloatSocket
        Amount to mix from A to B
    i.index : IntegerSocket
        `Index` on the geometries to sample from

    Outputs
    -------
    o.value : FloatSocket
        The final mixed float
    """

    _name = "Sample Mix Float"
    _asset_name = "Sample Mix Float"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "description": "Sample a float value from two different geometries and mix from A to B",
        "node_tool_idname": "geometry.sample_mix_float",
    }

    class _Inputs(SocketAccessor):
        a: GeometrySocket
        """Geometry A to sample and mix from"""
        b: GeometrySocket
        """Geometry B to sample and mix to"""
        value: FloatSocket
        """Float field to sample and mix"""
        factor: FloatSocket
        """Amount to mix from A to B"""
        index: IntegerSocket
        """`Index` on the geometries to sample from"""

    class _Outputs(SocketAccessor):
        value: FloatSocket
        """The final mixed float"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        a: InputGeometry = None,
        b: InputGeometry = None,
        value: InputFloat = 0.0,
        factor: InputFloat = 0.5,
        index: InputInteger = 0,
    ):
        super().__init__(
            **{"A": a, "B": b, "Value": value, "Factor": factor, "Index": index}
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        a = tree.inputs.geometry("A", description="Geometry A to sample and mix from")
        b = tree.inputs.geometry("B", description="Geometry B to sample and mix to")
        value = tree.inputs.float(
            "Value", 0.0, description="Float field to sample and mix", hide_value=True
        )
        factor = tree.inputs.float(
            "Factor",
            0.5,
            description="Amount to mix from A to B",
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
        value_1 = tree.outputs.float("Value", description="The final mixed float")

        sample_index = g.SampleIndex(geometry=a, value=value, index=index)
        sample_index_1 = g.SampleIndex(geometry=b, value=value, index=index)
        mix = g.Mix(
            factor_float=factor,
            a_float=sample_index,
            b_float=sample_index_1,
            clamp_factor=True,
        )
        vector = g.Vector(vector=(0.0, 0.0, 1.0))
        mix_1 = g.Mix(
            factor_float=factor,
            a_rotation=g.AxisAngleToRotation(axis=vector, angle=sample_index),
            b_rotation=g.AxisAngleToRotation(axis=vector, angle=sample_index_1),
            data_type="ROTATION",
            clamp_factor=True,
        )
        _rotation_to_axis_angle = g.RotationToAxisAngle(
            rotation=mix_1.o.result_rotation
        )

        mix >> value_1


ASSET = SampleMixFloat

ASSET_METADATA = {
    "description": "Sample a float value from two different geometries and mix from A to B",
    "catalog_id": "dd5f0199-fa8b-4b01-a972-2dc586a3e60f",
}
