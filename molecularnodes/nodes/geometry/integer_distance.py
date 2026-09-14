# Node-group asset "Integer Distance" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputInteger


class IntegerDistance(AssetGeometryGroup):
    """
    Integer Distance

    Parameters
    ----------
    a : InputInteger
        A
    b : InputInteger
        B
    distance : InputInteger
        Distance

    Inputs
    ------
    i.a : IntegerSocket
        A
    i.b : IntegerSocket
        B
    i.distance : IntegerSocket
        Distance

    Outputs
    -------
    o.cutoff : BooleanSocket
        Cutoff
    o.distance : IntegerSocket
        Distance
    """

    _name = "Integer Distance"
    _asset_name = "Integer Distance"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"

    class _Inputs(SocketAccessor):
        a: IntegerSocket
        """A"""
        b: IntegerSocket
        """B"""
        distance: IntegerSocket
        """Distance"""

    class _Outputs(SocketAccessor):
        cutoff: BooleanSocket
        """Cutoff"""
        distance: IntegerSocket
        """Distance"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        a: InputInteger = 0,
        b: InputInteger = 0,
        distance: InputInteger = 2,
    ):
        super().__init__(**{"A": a, "B": b, "Distance": distance})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        a = tree.inputs.integer("A", 0)
        b = tree.inputs.integer("B", 0)
        distance = tree.inputs.integer("Distance", 2)
        cutoff = tree.outputs.boolean("Cutoff")
        distance_1 = tree.outputs.integer("Distance")

        integer_math = abs(a - b)
        (integer_math >= distance) >> cutoff

        integer_math >> distance_1


ASSET = IntegerDistance

ASSET_METADATA = {
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
