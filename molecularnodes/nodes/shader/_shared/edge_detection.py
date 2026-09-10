# Node group 'Edge Detection' (ShaderNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
import math
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy import shader as s
from nodebpy.builder import CustomShaderGroup, FloatSocket, SocketAccessor
from nodebpy.types import InputFloat
from .offset_raycast import OffsetRaycast


class EdgeDetection(CustomShaderGroup):
    """
    Edge Detection

    Parameters
    ----------
    offset : InputFloat
        Offset

    Inputs
    ------
    i.offset : FloatSocket
        Offset

    Outputs
    -------
    o.co_planar_delta : FloatSocket
        Co-Planar Delta
    o.normal_delta : FloatSocket
        Normal Delta
    o.object_edge : FloatSocket
        Object Edge
    """

    _name = "Edge Detection"

    class _Inputs(SocketAccessor):
        offset: FloatSocket
        """Offset"""

    class _Outputs(SocketAccessor):
        co_planar_delta: FloatSocket
        """Co-Planar Delta"""
        normal_delta: FloatSocket
        """Normal Delta"""
        object_edge: FloatSocket
        """Object Edge"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        offset: InputFloat = 1.0,
    ):
        super().__init__(**{"Offset": offset})

    def _build_group(self, tree):
        offset = tree.inputs.float("Offset", 1.0, min_value=0.0, max_value=1.0)
        co_planar_delta = tree.outputs.float("Co-Planar Delta")
        normal_delta = tree.outputs.float("Normal Delta")
        object_edge = tree.outputs.float("Object Edge")

        geometry = s.Geometry()
        geometry_1 = s.Geometry()
        repeat_zone = g.RepeatZone(8)
        max_distance = repeat_zone.items.float("Max Distance")
        max_normal_delta = repeat_zone.items.float("Max Normal Delta")
        object_edge_1 = repeat_zone.items.float("Object Edge")
        map_range = g.WhiteNoiseTexture(
            vector=s.Geometry().o.position, noise_dimensions="4D"
        ).o.value.map_range(to_max=360.0)
        vector_rotate = g.VectorRotate.z_axis(
            g.CombineXYZ(x=1.0), angle=repeat_zone.iteration * 45.0 + map_range
        )
        vector = vector_rotate.o.vector * offset
        group = OffsetRaycast(x_offset=vector.x, y_offset=vector.y, length=10.0)
        math_1 = g.Math.greater_than(group.o.hit_distance, 0.0)
        _math_2 = 1.0 - group.o.self_hit
        vector_math = geometry.o.normal.dot(
            geometry.o.position + group.o.ray_direction * group.o.hit_distance
        )
        vector_math_1 = group.o.hit_normal.dot(
            geometry_1.o.normal * geometry_1.o.backfacing.mix.float(1.0, -1.0)
        )
        mix = g.Mix(
            factor_float=math_1,
            b_float=vector_math_1.acos().map_range(from_max=math.pi),
            clamp_factor=True,
        )
        (
            (geometry.o.position.dot(geometry.o.normal) - vector_math).max(
                max_distance.current
            )
            >> max_distance.next
        )
        max_normal_delta.current.max(mix.o.result_float) >> max_normal_delta.next
        (
            g.Mix(factor_float=math_1, clamp_factor=True).o.result_float.max(
                object_edge_1.current
            )
            >> object_edge_1.next
        )

        max_distance.result >> co_planar_delta
        max_normal_delta.result >> normal_delta
        object_edge_1.result >> object_edge
