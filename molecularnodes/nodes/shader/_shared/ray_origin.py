# Node group "Ray Origin" (ShaderNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import ShaderNodeTree
from nodebpy import TreeBuilder
from nodebpy import shader as s
from nodebpy.builder import CustomShaderGroup, SocketAccessor, VectorSocket


class RayOrigin(CustomShaderGroup):
    """
    Ray Origin

    Outputs
    -------
    o.vector : VectorSocket
        Vector
    """

    _name = "Ray Origin"

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        vector: VectorSocket
        """Vector"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self):
        super().__init__()

    def _build_group(self, tree: TreeBuilder[ShaderNodeTree]) -> None:
        vector = tree.outputs.vector("Vector")

        (
            s.Geometry().o.incoming * s.LightPath().o.ray_length
            + s.Geometry().o.position
            >> vector
        )
