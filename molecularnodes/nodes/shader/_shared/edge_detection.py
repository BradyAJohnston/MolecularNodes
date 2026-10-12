# Node group "Edge Detection" (ShaderNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
import math
from typing import TYPE_CHECKING
from bpy.types import ShaderNodeTree
from nodebpy import TreeBuilder
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

    def _build_group(self, tree: TreeBuilder[ShaderNodeTree]) -> None:
        offset = tree.inputs.float("Offset", 1.0, min_value=0.0, max_value=1.0)
        co_planar_delta = tree.outputs.float("Co-Planar Delta")
        normal_delta = tree.outputs.float("Normal Delta")
        object_edge = tree.outputs.float("Object Edge")

        geometry = s.Geometry()
        geometry_1 = s.Geometry()
        repeat_zone = s.RepeatZone(8)
        max_distance = repeat_zone.items.float(name="Max Distance")
        max_normal_delta = repeat_zone.items.float(name="Max Normal Delta")
        object_edge_1 = repeat_zone.items.float(name="Object Edge")
        vector_math = (
            geometry_1.o.normal
            * s.Mix.float(geometry_1.o.backfacing, 1.0, -1.0).o.result_float
        )
        map_range = s.WhiteNoiseTexture(
            vector=s.Geometry().o.position, noise_dimensions="4D"
        ).o.value.map_range(to_max=360.0)
        vector_math_1 = (
            s.VectorRotate.z_axis(
                s.CombineXYZ(1.0), angle=repeat_zone.iteration * 45.0 + map_range
            ).o.vector
            * offset
        )
        offset_raycast = OffsetRaycast(
            x_offset=vector_math_1.x, y_offset=vector_math_1.y, length=10.0
        )
        math_1 = s.Math.greater_than(offset_raycast.o.hit_distance, 0.0)
        _math_2 = 1.0 - offset_raycast.o.self_hit
        vector_math_2 = geometry.o.normal.dot(
            geometry.o.position
            + offset_raycast.o.ray_direction * offset_raycast.o.hit_distance
        )
        mix = math_1.o.value.mix.float(
            0.0,
            offset_raycast.o.hit_normal.dot(vector_math)
            .acos()
            .map_range(from_max=math.pi),
        )
        (
            (geometry.o.position.dot(geometry.o.normal) - vector_math_2).max(
                max_distance.current
            )
            >> max_distance.next
        )
        max_normal_delta.current.max(mix) >> max_normal_delta.next
        (
            math_1.o.value.mix.float(0.0, 0.0).max(object_edge_1.current)
            >> object_edge_1.next
        )

        max_distance.result >> co_planar_delta
        max_normal_delta.result >> normal_delta
        object_edge_1.result >> object_edge
