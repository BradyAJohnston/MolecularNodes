# Node group "Offset Raycast" (ShaderNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import ShaderNodeTree
from nodebpy import TreeBuilder
from nodebpy import shader as s
from nodebpy.builder import CustomShaderGroup, FloatSocket, SocketAccessor, VectorSocket
from nodebpy.types import InputFloat
from .ray_origin import RayOrigin
from .ray_tangent import RayTangent


class OffsetRaycast(CustomShaderGroup):
    """
    Offset Raycast

    Parameters
    ----------
    x_offset : InputFloat
        X Offset
    y_offset : InputFloat
        Y Offset
    length : InputFloat
        Length

    Inputs
    ------
    i.x_offset : FloatSocket
        X Offset
    i.y_offset : FloatSocket
        Y Offset
    i.length : FloatSocket
        Length

    Outputs
    -------
    o.is_hit : FloatSocket
        Is Hit
    o.self_hit : FloatSocket
        Self Hit
    o.hit_position : VectorSocket
        Hit Position
    o.hit_normal : VectorSocket
        Hit Normal
    o.hit_distance : FloatSocket
        Hit Distance
    o.ray_direction : VectorSocket
        Ray Direction
    """

    _name = "Offset Raycast"

    class _Inputs(SocketAccessor):
        x_offset: FloatSocket
        """X Offset"""
        y_offset: FloatSocket
        """Y Offset"""
        length: FloatSocket
        """Length"""

    class _Outputs(SocketAccessor):
        is_hit: FloatSocket
        """Is Hit"""
        self_hit: FloatSocket
        """Self Hit"""
        hit_position: VectorSocket
        """Hit Position"""
        hit_normal: VectorSocket
        """Hit Normal"""
        hit_distance: FloatSocket
        """Hit Distance"""
        ray_direction: VectorSocket
        """Ray Direction"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        x_offset: InputFloat = 1.0,
        y_offset: InputFloat = 1.0,
        length: InputFloat = 0.0,
    ):
        super().__init__(
            **{"X Offset": x_offset, "Y Offset": y_offset, "Length": length}
        )

    def _build_group(self, tree: TreeBuilder[ShaderNodeTree]) -> None:
        x_offset = tree.inputs.float(
            "X Offset", 1.0, min_value=-10_000.0, max_value=10_000.0
        )
        y_offset = tree.inputs.float(
            "Y Offset", 1.0, min_value=-10_000.0, max_value=10_000.0
        )
        length = tree.inputs.float("Length", 0.0)
        is_hit = tree.outputs.float("Is Hit")
        self_hit = tree.outputs.float("Self Hit")
        hit_position = tree.outputs.vector("Hit Position")
        hit_normal = tree.outputs.vector("Hit Normal")
        hit_distance = tree.outputs.float("Hit Distance")
        ray_direction = tree.outputs.vector("Ray Direction")

        group = RayTangent()
        group_1 = RayOrigin()
        vector_math = (
            s.Geometry().o.position
            + (group.o.tangent * x_offset + group.o.bitangent * y_offset)
            - group_1
        )
        vector_math_1 = vector_math.normalize()
        vector_math_2 = vector_math * 0.5
        raycast = s.Raycast(
            position=group_1.o.vector + vector_math_2,
            direction=vector_math_1,
            length=s.LightPath().o.ray_length + length,
        )
        (
            raycast.o.hit_distance + vector_math_2.length() - s.LightPath().o.ray_length
            >> hit_distance
        )

        raycast >> is_hit
        raycast.o.self_hit >> self_hit
        raycast.o.hit_position >> hit_position
        raycast.o.hit_normal >> hit_normal
        vector_math_1 >> ray_direction
