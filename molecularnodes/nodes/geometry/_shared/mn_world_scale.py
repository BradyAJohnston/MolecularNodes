# Node group ".MN_world_scale" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import CustomGeometryGroup, FloatSocket, SocketAccessor


class MN_world_scale(CustomGeometryGroup):
    """
    .MN_world_scale

    Outputs
    -------
    o.world_scale : FloatSocket
        world_scale
    """

    _name = ".MN_world_scale"
    _tree_properties = {"node_tool_idname": "geometry._mn_world_scale"}

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        world_scale: FloatSocket

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self):
        super().__init__()

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        world_scale = tree.outputs.float("world_scale", 0.01)

        g.Value(0.1) >> world_scale
