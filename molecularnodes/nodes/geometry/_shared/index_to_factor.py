# Node group 'Index to Factor' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    CustomGeometryGroup,
    FloatSocket,
    IntegerSocket,
    SocketAccessor,
)
from nodebpy.types import InputInteger


class IndexToFactor(CustomGeometryGroup):
    """
    Index to Factor

    Parameters
    ----------
    index : InputInteger
        The index starting at 0 within the overall group of size `Size`
    size : InputInteger
        The size of the group that will be used to compute the factor.

    Inputs
    ------
    i.index : IntegerSocket
        The index starting at 0 within the overall group of size `Size`
    i.size : IntegerSocket
        The size of the group that will be used to compute the factor.

    Outputs
    -------
    o.factor : FloatSocket
        Factor
    """

    _name = "Index to Factor"
    _color_tag = "CONVERTER"

    class _Inputs(SocketAccessor):
        index: IntegerSocket
        """The index starting at 0 within the overall group of size `Size`"""
        size: IntegerSocket
        """The size of the group that will be used to compute the factor."""

    class _Outputs(SocketAccessor):
        factor: FloatSocket
        """Factor"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        index: InputInteger = 0,
        size: InputInteger = 0,
    ):
        super().__init__(**{"Index": index, "Size": size})

    def _build_group(self, tree):
        index = tree.inputs.integer(
            "Index",
            0,
            description="The index starting at 0 within the overall group of size `Size`",
        )
        size = tree.inputs.integer(
            "Size",
            0,
            description="The size of the group that will be used to compute the factor.",
        )
        factor = tree.outputs.float("Factor", subtype="FACTOR")

        map_range = g.MapRange(value=index, from_max=size - 1, clamp=True)

        map_range >> factor
