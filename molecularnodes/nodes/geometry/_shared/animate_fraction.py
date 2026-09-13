# Node group "Animate Fraction" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    BooleanSocket,
    CustomGeometryGroup,
    FloatSocket,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputFloat


class AnimateFraction(CustomGeometryGroup):
    """
    Interpolate the fraction component of a float

    Parameters
    ----------
    interpolate : InputBoolean
        Interpolate
    smoother_step : InputBoolean
        Smoother Step
    float : InputFloat
        Float

    Inputs
    ------
    i.interpolate : BooleanSocket
        Interpolate
    i.smoother_step : BooleanSocket
        Smoother Step
    i.float : FloatSocket
        Float

    Outputs
    -------
    o.float : FloatSocket
        Float
    """

    _name = "Animate Fraction"
    _color_tag = "CONVERTER"
    _tree_properties = {
        "description": "Interpolate the fraction component of a float",
        "node_tool_idname": "geometry.animate_fraction",
    }

    class _Inputs(SocketAccessor):
        interpolate: BooleanSocket
        """Interpolate"""
        smoother_step: BooleanSocket
        """Smoother Step"""
        float: FloatSocket
        """Float"""

    class _Outputs(SocketAccessor):
        float: FloatSocket
        """Float"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        interpolate: InputBoolean = False,
        smoother_step: InputBoolean = False,
        float: InputFloat = 0.0,
    ):
        super().__init__(
            **{
                "Interpolate": interpolate,
                "Smoother Step": smoother_step,
                "Float": float,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        interpolate = tree.inputs.boolean("Interpolate", False)
        smoother_step = tree.inputs.boolean("Smoother Step", False)
        float = tree.inputs.float("Float", 0.0)
        float_1 = tree.outputs.float("Float")

        math_1 = float.fraction()
        (
            interpolate.switch.float(
                float.floor(),
                smoother_step.switch.float(
                    math_1, math_1.map_range(interpolation_type="SMOOTHERSTEP")
                ),
            )
            >> float_1
        )
