# Node group 'CA Value Vector' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import CustomGeometryGroup, SocketAccessor, VectorSocket
from nodebpy.types import InputVector
from ..atom_name import AtomName


class CAValueVector(CustomGeometryGroup):
    """
    CA Value Vector

    Parameters
    ----------
    vector : InputVector
        Vector

    Inputs
    ------
    i.vector : VectorSocket
        Vector

    Outputs
    -------
    o.vector : VectorSocket
        Vector
    """

    _name = "CA Value Vector"
    _color_tag = "CONVERTER"

    class _Inputs(SocketAccessor):
        vector: VectorSocket
        """Vector"""

    class _Outputs(SocketAccessor):
        vector: VectorSocket
        """Vector"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        vector: InputVector = None,
    ):
        super().__init__(**{"Vector": vector})

    def _build_group(self, tree):
        vector = tree.inputs.vector("Vector", (0.0, 0.0, 0.0), hide_value=True)
        vector_1 = tree.outputs.vector("Vector")

        (
            g.IndexSwitch.vector(AtomName(), ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), vector))
            >> vector_1
        )
