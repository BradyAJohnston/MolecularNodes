# Node group "Lerp Position" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    CustomGeometryGroup,
    FloatSocket,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputFloat, InputVector


class LerpPosition(CustomGeometryGroup):
    """
    Lerp Position

    Parameters
    ----------
    a : InputVector
        a
    b : InputVector
        b
    deltat : InputFloat
        deltaT
    decay : InputFloat
        Decay

    Inputs
    ------
    i.a : VectorSocket
        a
    i.b : VectorSocket
        b
    i.deltat : FloatSocket
        deltaT
    i.decay : FloatSocket
        Decay

    Outputs
    -------
    o.position : VectorSocket
        Position
    o.length : FloatSocket
        Length
    """

    _name = "Lerp Position"
    _color_tag = "CONVERTER"

    class _Inputs(SocketAccessor):
        a: VectorSocket
        b: VectorSocket
        deltat: FloatSocket
        """deltaT"""
        decay: FloatSocket
        """Decay"""

    class _Outputs(SocketAccessor):
        position: VectorSocket
        """Position"""
        length: FloatSocket
        """Length"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        a: InputVector = None,
        b: InputVector = None,
        deltat: InputFloat = 0.5,
        decay: InputFloat = 0.5,
    ):
        super().__init__(**{"a": a, "b": b, "deltaT": deltat, "Decay": decay})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        a = tree.inputs.vector(
            "a",
            (0.0, 0.0, 0.0),
            min_value=-10_000.0,
            max_value=10_000.0,
            default_input="POSITION",
        )
        b = tree.inputs.vector(
            "b", (0.0, 0.0, 0.0), min_value=-10_000.0, max_value=10_000.0
        )
        deltat = tree.inputs.float(
            "deltaT",
            0.5,
            min_value=-10_000.0,
            max_value=10_000.0,
            structure_type="SINGLE",
            force_non_field=True,
        )
        decay = tree.inputs.float("Decay", 0.5, min_value=-10_000.0, max_value=10_000.0)
        position = tree.outputs.vector("Position")
        length = tree.outputs.float("Length")

        vector_math = a - b
        b + vector_math * g.Math.exponent(deltat * (decay * -1.0)) >> position
        vector_math.length() >> length
