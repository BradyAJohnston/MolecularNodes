# Node group ".Check End Face Corner" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
    IntegerSocket,
    SocketAccessor,
)
from nodebpy.types import InputInteger


class CheckEndFaceCorner(CustomGeometryGroup):
    """
    .Check End Face Corner

    Parameters
    ----------
    captured_index : InputInteger
        Captured Index

    Inputs
    ------
    i.captured_index : IntegerSocket
        Captured Index

    Outputs
    -------
    o.is_end_face_corner : BooleanSocket
        Is End Face Corner
    """

    _name = ".Check End Face Corner"
    _color_tag = "CONVERTER"

    class _Inputs(SocketAccessor):
        captured_index: IntegerSocket
        """Captured Index"""

    class _Outputs(SocketAccessor):
        is_end_face_corner: BooleanSocket
        """Is End Face Corner"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        captured_index: InputInteger = 0,
    ):
        super().__init__(**{"Captured Index": captured_index})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        captured_index = tree.inputs.integer("Captured Index", 0)
        is_end_face_corner = tree.outputs.boolean("Is End Face Corner")

        (
            (
                g.Compare.integer.equal(captured_index, 0).o.result
                & g.Compare.integer.not_equal(
                    captured_index.corner.at(g.OffsetCornerInFace(offset=2)), 1
                )
            )
            >> is_end_face_corner
        )
