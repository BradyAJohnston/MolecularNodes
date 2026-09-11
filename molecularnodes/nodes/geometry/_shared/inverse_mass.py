# Node group "Inverse Mass" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    BooleanSocket,
    CustomGeometryGroup,
    FloatSocket,
    SocketAccessor,
)
from ..fallback_float import FallbackFloat


class InverseMass(CustomGeometryGroup):
    """
    Inverse Mass

    Outputs
    -------
    o.exists : BooleanSocket
        Exists
    o.w : FloatSocket
        w
    o.mass : FloatSocket
        mass
    """

    _name = "Inverse Mass"
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        exists: BooleanSocket
        """Exists"""
        w: FloatSocket
        mass: FloatSocket

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self):
        super().__init__()

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        exists = tree.outputs.boolean("Exists")
        w = tree.outputs.float("w")
        mass = tree.outputs.float("mass")

        group = FallbackFloat(name="mass", fallback=1.0)
        string = g.String(string="inverse_mass")
        FallbackFloat(name=string, fallback=1.0 / group) >> w
        named_attribute = g.NamedAttribute.float(string)

        named_attribute.o.exists >> exists
        group >> mass
