# Node-group asset "Animate Ease" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
import math
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputFloat, InputMenu
from ._shared.ease_switch import EaseSwitch


class AnimateEase(AssetGeometryGroup):
    """
    Animate Ease

    Parameters
    ----------
    value : InputFloat
        Linear progress to ease, from 0 to 1. Link `Animate Value` or `Animate Stagger` to animate it
    interpolation : InputMenu | Literal["Linear", "Sinusoidal", "Quadratic", "Cubic", "Quartic", "Quintic", "Exponential", "Circular", "Back", "Bounce", "Elastic"]
        Shape of the easing curve, following Robert Penner's easing functions (easings.net)
    ease : InputMenu | Literal["In", "Out", "In Out"]
        Apply the curve at the start (In), the end (Out) or both ends (In Out) of the transition
    clamp : InputBoolean
        Clamp the input value to 0..1 before easing so the output never leaves the From..To range
    from_ : InputFloat
        Output value when the input is 0
    to : InputFloat
        Output value when the input is 1
    overshoot : InputFloat
        How far the `Back` curve overshoots before settling
    period : InputFloat
        Period of the `Elastic` oscillation as a fraction of the transition
    amplitude : InputFloat
        Amplitude of the `Elastic` oscillation. Values below 1 behave as 1

    Inputs
    ------
    i.value : FloatSocket
        Linear progress to ease, from 0 to 1. Link `Animate Value` or `Animate Stagger` to animate it
    i.interpolation : MenuSocket
        Shape of the easing curve, following Robert Penner's easing functions (easings.net)
    i.ease : MenuSocket
        Apply the curve at the start (In), the end (Out) or both ends (In Out) of the transition
    i.clamp : BooleanSocket
        Clamp the input value to 0..1 before easing so the output never leaves the From..To range
    i.from_ : FloatSocket
        Output value when the input is 0
    i.to : FloatSocket
        Output value when the input is 1
    i.overshoot : FloatSocket
        How far the `Back` curve overshoots before settling
    i.period : FloatSocket
        Period of the `Elastic` oscillation as a fraction of the transition
    i.amplitude : FloatSocket
        Amplitude of the `Elastic` oscillation. Values below 1 behave as 1

    Outputs
    -------
    o.value : FloatSocket
        Eased value between `From` and `To`
    """

    _name = "Animate Ease"
    _asset_name = "Animate Ease"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"

    class _Inputs(SocketAccessor):
        value: FloatSocket
        """Linear progress to ease, from 0 to 1. Link `Animate Value` or `Animate Stagger` to animate it"""
        interpolation: MenuSocket
        """Shape of the easing curve, following Robert Penner's easing functions (easings.net)"""
        ease: MenuSocket
        """Apply the curve at the start (In), the end (Out) or both ends (In Out) of the transition"""
        clamp: BooleanSocket
        """Clamp the input value to 0..1 before easing so the output never leaves the From..To range"""
        from_: FloatSocket
        """Output value when the input is 0"""
        to: FloatSocket
        """Output value when the input is 1"""
        overshoot: FloatSocket
        """How far the `Back` curve overshoots before settling"""
        period: FloatSocket
        """Period of the `Elastic` oscillation as a fraction of the transition"""
        amplitude: FloatSocket
        """Amplitude of the `Elastic` oscillation. Values below 1 behave as 1"""

    class _Outputs(SocketAccessor):
        value: FloatSocket
        """Eased value between `From` and `To`"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        value: InputFloat = 0.0,
        interpolation: InputMenu
        | Literal[
            "Linear",
            "Sinusoidal",
            "Quadratic",
            "Cubic",
            "Quartic",
            "Quintic",
            "Exponential",
            "Circular",
            "Back",
            "Bounce",
            "Elastic",
        ] = "Cubic",
        ease: InputMenu | Literal["In", "Out", "In Out"] = "In Out",
        clamp: InputBoolean = True,
        from_: InputFloat = 0.0,
        to: InputFloat = 1.0,
        overshoot: InputFloat = 1.70158,
        period: InputFloat = 0.3,
        amplitude: InputFloat = 1.0,
    ):
        super().__init__(
            **{
                "Value": value,
                "Interpolation": interpolation,
                "Ease": ease,
                "Clamp": clamp,
                "From": from_,
                "To": to,
                "Overshoot": overshoot,
                "Period": period,
                "Amplitude": amplitude,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        value = tree.inputs.float(
            "Value",
            0.0,
            description="Linear progress to ease, from 0 to 1. Link `Animate Value` or `Animate Stagger` to animate it",
        )
        interpolation = tree.inputs.menu(
            "Interpolation",
            description="Shape of the easing curve, following Robert Penner's easing functions (easings.net)",
        )
        ease = tree.inputs.menu(
            "Ease",
            description="Apply the curve at the start (In), the end (Out) or both ends (In Out) of the transition",
        )
        clamp = tree.inputs.boolean(
            "Clamp",
            True,
            description="Clamp the input value to 0..1 before easing so the output never leaves the From..To range",
        )
        from_ = tree.inputs.float(
            "From",
            0.0,
            description="Output value when the input is 0",
            min_value=-10_000.0,
            max_value=10_000.0,
        )
        to = tree.inputs.float(
            "To",
            1.0,
            description="Output value when the input is 1",
            min_value=-10_000.0,
            max_value=10_000.0,
        )
        with tree.inputs.panel("Curve", default_closed=True):
            overshoot = tree.inputs.float(
                "Overshoot",
                1.70158,
                description="How far the `Back` curve overshoots before settling",
                min_value=0.0,
                max_value=10.0,
            )
            period = tree.inputs.float(
                "Period",
                0.3,
                description="Period of the `Elastic` oscillation as a fraction of the transition",
                min_value=0.01,
                max_value=10.0,
            )
            amplitude = tree.inputs.float(
                "Amplitude",
                1.0,
                description="Amplitude of the `Elastic` oscillation. Values below 1 behave as 1",
                min_value=0.0,
                max_value=10.0,
            )
        value_1 = tree.outputs.float(
            "Value", description="Eased value between `From` and `To`"
        )

        switch = clamp.switch.float(value, value.clamp())
        with g.Frame("Ease argument"):
            ease_switch = EaseSwitch(
                ease=ease, in_=overshoot, out=overshoot, in_out=overshoot * 1.525
            )
            ease_switch_1 = EaseSwitch(
                ease=ease, in_=period, out=period, in_out=period * 1.5
            )
            _string = g.String(
                string="Each curve is built once in its ease-in form f(x). Out is 1 - f(1 - t) and In Out evaluates f(2t) / 2 for the first half and 1 - f(2 - 2t) / 2 for the second, so this block picks the argument x and the block after the curves maps f(x) back. Back and Elastic use Penner's larger overshoot and period for In Out. The Ease menu goes through the Ease Switch group because one menu socket feeding several Menu Switch nodes loses its items on a headless build."
            )
            compare = switch > 0.5
            ease_switch_2 = EaseSwitch(
                ease=ease,
                in_=switch,
                out=1.0 - switch,
                in_out=compare.switch.float(switch * 2.0, 2.0 - switch * 2.0),
            )
        with g.Frame("Ease-in curves"):
            _string_1 = g.String(
                string="Reference: easings.net. x is the ease argument in 0..1. Exponential and Elastic are pinned to exactly 0 at x = 0 (and Elastic to 1 at x = 1) as in the reference, since 2^-10 is not zero."
            )
            math_1 = 1.0 - (ease_switch_2.o.value * (math.pi / 2)).cos()
            math_2 = ease_switch_2.o.value * ease_switch_2
            math_3 = ease_switch_2.o.value**3.0
            math_4 = ease_switch_2.o.value**4.0
            math_5 = ease_switch_2.o.value**5.0
            math_6 = 1.0 - (1.0 - ease_switch_2.o.value * ease_switch_2).sqrt()
            math_7 = (
                (ease_switch.o.value + 1.0) * ease_switch_2.o.value** 3.0
                - ease_switch.o.value * ease_switch_2 * ease_switch_2
            )
            switch_1 = (ease_switch_2 <= 0.0).switch.float(
                2.0 ** (ease_switch_2.o.value * 10.0 - 10.0)
            )
        with g.Frame("Elastic"):
            math_8 = amplitude.max(1.0)
            _string_2 = g.String(
                string="Penner's elastic with amplitude a (at least 1) and period p: s = p / (2 pi) * asin(1 / a), f(x) = -a * 2^(10x - 10) * sin((x - 1 - s) * 2 pi / p)."
            )
            math_9 = (
                ease_switch_2.o.value
                - 1.0
                - (1.0 / math_8).asin() * (ease_switch_1.o.value / math.tau)
            ) / ease_switch_1
            switch_2 = (ease_switch_2 >= 1.0).switch.float(
                math_8
                * 2.0 ** (ease_switch_2.o.value * 10.0 - 10.0)
                * (math_9 * math.tau).sin()
                * -1.0,
                1.0,
            )
            switch_3 = (ease_switch_2 <= 0.0).switch.float(switch_2)
        with g.Frame("Bounce"):
            _string_3 = g.String(
                string="Penner defines Bounce by its ease-out form b(u) with four parabolic arcs; the ease-in form is 1 - b(1 - x)."
            )
            math_10 = 1.0 - ease_switch_2
            math_11 = math_10 - 0.5454546
            math_12 = math_10 - 0.8181818
            math_13 = math_10 - 0.9545454
            switch_4 = (math_10 < 0.909091).switch.float(
                math_13 * math_13 * 7.5625 + 0.984375,
                math_12 * math_12 * 7.5625 + 0.9375,
            )
            switch_5 = (math_10 < 0.3636364).switch.float(
                (math_10 < 0.7272727).switch.float(
                    switch_4, math_11 * math_11 * 7.5625 + 0.75
                ),
                math_10 * math_10 * 7.5625,
            )
            math_14 = 1.0 - switch_5
        menu_switch = g.MenuSwitch.float(
            interpolation,
            {
                "Linear": ease_switch_2,
                "Sinusoidal": math_1,
                "Quadratic": math_2,
                "Cubic": math_3,
                "Quartic": math_4,
                "Quintic": math_5,
                "Exponential": switch_1,
                "Circular": math_6,
                "Back": math_7,
                "Bounce": math_14,
                "Elastic": switch_3,
            },
        )
        with g.Frame("Ease direction"):
            ease_switch_3 = EaseSwitch(
                ease=ease,
                in_=menu_switch.o.output,
                out=1.0 - menu_switch.o.output,
                in_out=compare.switch.float(
                    menu_switch.o.output * 0.5, 1.0 - menu_switch.o.output * 0.5
                ),
            )
        from_ + (to - from_) * ease_switch_3 >> value_1

        interpolation.default_value = "Cubic"
        ease.default_value = "In Out"


ASSET = AnimateEase

ASSET_METADATA = {
    "description": "Ease a 0..1 value with Penner easing curves (In, Out or In Out) and map it to a From..To range",
    "catalog_id": "85730213-4c2e-469f-b333-52ac53adf274",
}
