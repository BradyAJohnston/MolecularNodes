# Node group "Ray Tangent" (ShaderNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import ShaderNodeTree
from nodebpy import TreeBuilder
from nodebpy import shader as s
from nodebpy.builder import CustomShaderGroup, SocketAccessor, VectorSocket


class RayTangent(CustomShaderGroup):
    """
    Ray Tangent

    Outputs
    -------
    o.tangent : VectorSocket
        Tangent
    o.bitangent : VectorSocket
        Bitangent
    """

    _name = "Ray Tangent"

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        tangent: VectorSocket
        """Tangent"""
        bitangent: VectorSocket
        """Bitangent"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self):
        super().__init__()

    def _build_group(self, tree: TreeBuilder[ShaderNodeTree]) -> None:
        tangent = tree.outputs.vector("Tangent")
        bitangent = tree.outputs.vector("Bitangent")

        geometry = s.Geometry()
        vector_transform = s.VectorTransform(
            vector=(0.0, 1.0, 0.0),
            vector_type="NORMAL",
            convert_from="CAMERA",
            convert_to="WORLD",
        )
        vector_math = (geometry.o.incoming * -1.0).cross(vector_transform).normalize()
        geometry.o.incoming.cross(vector_math) >> bitangent

        vector_math >> tangent
