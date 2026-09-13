# Node group "Constraint Distance" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    CustomGeometryGroup,
    FloatSocket,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputFloat, InputVector


class ConstraintDistance(CustomGeometryGroup):
    """
    Constraint Distance

    Parameters
    ----------
    target : InputVector
        Target
    self_2 : InputVector
        Self
    distance : InputFloat
        Distance
    w1 : InputFloat
        W1
    w2 : InputFloat
        W2
    alpha : InputFloat
        alpha
    deltat : InputFloat
        deltaT

    Inputs
    ------
    i.target : VectorSocket
        Target
    i.self_2 : VectorSocket
        Self
    i.distance : FloatSocket
        Distance
    i.w1 : FloatSocket
        W1
    i.w2 : FloatSocket
        W2
    i.alpha : FloatSocket
        alpha
    i.deltat : FloatSocket
        deltaT

    Outputs
    -------
    o.correction : VectorSocket
        Correction
    o.value : FloatSocket
        Value
    """

    _name = "Constraint Distance"
    _color_tag = "VECTOR"

    class _Inputs(SocketAccessor):
        target: VectorSocket
        """Target"""
        self_2: VectorSocket
        """Self"""
        distance: FloatSocket
        """Distance"""
        w1: FloatSocket
        """W1"""
        w2: FloatSocket
        """W2"""
        alpha: FloatSocket
        deltat: FloatSocket
        """deltaT"""

    class _Outputs(SocketAccessor):
        correction: VectorSocket
        """Correction"""
        value: FloatSocket
        """Value"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        target: InputVector = None,
        self_2: InputVector = None,
        distance: InputFloat = 0.5,
        w1: InputFloat = 0.5,
        w2: InputFloat = 0.5,
        alpha: InputFloat = 0.0,
        deltat: InputFloat = 0.0,
    ):
        super().__init__(
            **{
                "Target": target,
                "Self": self_2,
                "Distance": distance,
                "W1": w1,
                "W2": w2,
                "alpha": alpha,
                "deltaT": deltat,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        target = tree.inputs.vector(
            "Target", (0.0, 0.0, 0.0), min_value=-10_000.0, max_value=10_000.0
        )
        self = tree.inputs.vector(
            "Self",
            (0.0, 0.0, 0.0),
            min_value=-10_000.0,
            max_value=10_000.0,
            default_input="POSITION",
        )
        distance = tree.inputs.float(
            "Distance", 0.5, min_value=-10_000.0, max_value=10_000.0
        )
        w1 = tree.inputs.float("W1", 0.5, min_value=0.0, max_value=10_000.0)
        w2 = tree.inputs.float("W2", 0.5, min_value=0.0, max_value=10_000.0)
        alpha = tree.inputs.float("alpha", 0.0, min_value=-10_000.0, max_value=10_000.0)
        deltat = tree.inputs.float(
            "deltaT",
            0.0,
            min_value=-10_000.0,
            max_value=10_000.0,
            structure_type="SINGLE",
            force_non_field=True,
        )
        correction = tree.outputs.vector("Correction")
        value = tree.outputs.float("Value")

        vector_math = target - self
        vector_math_1 = vector_math.length()
        (
            vector_math.normalize()
            * (
                (vector_math_1 - distance)
                * (w1 / (w1 + w2 + alpha / (deltat * deltat)))
            )
            >> correction
        )

        vector_math_1 >> value
