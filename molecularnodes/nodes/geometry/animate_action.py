# Node-group asset "Animate Action" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    IntegerSocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import Default, InputBoolean, InputFloat, InputInteger, InputMenu
from .between_float import BetweenFloat
from .ease_value import EaseValue
from .stagger_value import StaggerValue


class AnimateAction(AssetGeometryGroup):
    """
    Animate Action

    Parameters
    ----------
    start : InputFloat
        Start
    length : InputFloat
        Length
    time : InputMenu | Literal["Seconds", "Frames", "Value"]
        Which source to use as time for animating the action
    value : InputFloat
        Becomes the output value if it is chosen by the menu input
    stagger : InputBoolean
        Stagger
    width : InputFloat
        How many `ID`s are staggered at the same time during the animation.
    id : InputInteger
        The `ID` over which we stagger the animation. When unconnected: The "id" attribute if available, otherwise the index.
    group_id : InputInteger
        Compute the staggered animation individually for each `Group ID`
    interpolation : InputMenu | Literal["Linear", "Sinusoidal", "Quadratic", "Cubic", "Quartic", "Quintic", "Exponential", "Circular", "Back", "Bounce", "Elastic"]
        Shape of the easing curve, following Robert Penner's easing functions (easings.net)
    ease : InputMenu | Literal["In", "Out", "In Out"]
        Apply the curve at the start (In), the end (Out) or both ends (In Out) of the transition
    overshoot : InputFloat
        How far the `Back` curve overshoots before settling
    period : InputFloat
        Period of the `Elastic` oscillation as a fraction of the transition
    amplitude : InputFloat
        Amplitude of the `Elastic` oscillation. Values below 1 behave as 1

    Inputs
    ------
    i.start : FloatSocket
        Start
    i.length : FloatSocket
        Length
    i.time : MenuSocket
        Which source to use as time for animating the action
    i.value : FloatSocket
        Becomes the output value if it is chosen by the menu input
    i.stagger : BooleanSocket
        Stagger
    i.width : FloatSocket
        How many `ID`s are staggered at the same time during the animation.
    i.id : IntegerSocket
        The `ID` over which we stagger the animation.
    i.group_id : IntegerSocket
        Compute the staggered animation individually for each `Group ID`
    i.interpolation : MenuSocket
        Shape of the easing curve, following Robert Penner's easing functions (easings.net)
    i.ease : MenuSocket
        Apply the curve at the start (In), the end (Out) or both ends (In Out) of the transition
    i.overshoot : FloatSocket
        How far the `Back` curve overshoots before settling
    i.period : FloatSocket
        Period of the `Elastic` oscillation as a fraction of the transition
    i.amplitude : FloatSocket
        Amplitude of the `Elastic` oscillation. Values below 1 behave as 1

    Outputs
    -------
    o.value : FloatSocket
        Value
    o.active : BooleanSocket
        Whether the input `Value` is between (and including) the lower and upper bounds
    o.end : FloatSocket
        End
    """

    _name = "Animate Action"
    _asset_name = "Animate Action"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        start: FloatSocket
        """Start"""
        length: FloatSocket
        """Length"""
        time: MenuSocket
        """Which source to use as time for animating the action"""
        value: FloatSocket
        """Becomes the output value if it is chosen by the menu input"""
        stagger: BooleanSocket
        """Stagger"""
        width: FloatSocket
        """How many `ID`s are staggered at the same time during the animation."""
        id: IntegerSocket
        """The `ID` over which we stagger the animation."""
        group_id: IntegerSocket
        """Compute the staggered animation individually for each `Group ID`"""
        interpolation: MenuSocket
        """Shape of the easing curve, following Robert Penner's easing functions (easings.net)"""
        ease: MenuSocket
        """Apply the curve at the start (In), the end (Out) or both ends (In Out) of the transition"""
        overshoot: FloatSocket
        """How far the `Back` curve overshoots before settling"""
        period: FloatSocket
        """Period of the `Elastic` oscillation as a fraction of the transition"""
        amplitude: FloatSocket
        """Amplitude of the `Elastic` oscillation. Values below 1 behave as 1"""

    class _Outputs(SocketAccessor):
        value: FloatSocket
        """Value"""
        active: BooleanSocket
        """Whether the input `Value` is between (and including) the lower and upper bounds"""
        end: FloatSocket
        """End"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        start: InputFloat = 1.0,
        length: InputFloat = 5.0,
        time: InputMenu | Literal["Seconds", "Frames", "Value"] = "Seconds",
        value: InputFloat = 0.0,
        stagger: InputBoolean = False,
        width: InputFloat = 10.0,
        id: InputInteger = Default.ID_OR_INDEX,
        group_id: InputInteger = 0,
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
        ] = "Linear",
        ease: InputMenu | Literal["In", "Out", "In Out"] = "In",
        overshoot: InputFloat = 1.70158,
        period: InputFloat = 0.3,
        amplitude: InputFloat = 1.0,
    ):
        super().__init__(
            **{
                "Start": start,
                "Length": length,
                "Time": time,
                "Value": value,
                "Stagger": stagger,
                "Width": width,
                "ID": id,
                "Group ID": group_id,
                "Interpolation": interpolation,
                "Ease": ease,
                "Overshoot": overshoot,
                "Period": period,
                "Amplitude": amplitude,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        start = tree.inputs.float("Start", 1.0, min_value=-10_000.0, max_value=10_000.0)
        length = tree.inputs.float("Length", 5.0, min_value=0.0)
        time = tree.inputs.menu(
            "Time",
            description="Which source to use as time for animating the action",
            optional_label=True,
        )
        value = tree.inputs.float(
            "Value",
            0.0,
            description="Becomes the output value if it is chosen by the menu input",
        )
        with tree.inputs.panel("Stagger", default_closed=True):
            stagger = tree.inputs.boolean("Stagger", False, is_panel_toggle=True)
            width = tree.inputs.float(
                "Width",
                10.0,
                description="How many `ID`s are staggered at the same time during the animation.",
                min_value=0.0,
                max_value=10_000.0,
            )
            id = tree.inputs.integer(
                "ID",
                0,
                description="The `ID` over which we stagger the animation.",
                hide_value=True,
                default_input="ID_OR_INDEX",
            )
            group_id = tree.inputs.integer(
                "Group ID",
                0,
                description="Compute the staggered animation individually for each `Group ID`",
                hide_value=True,
            )
        with tree.inputs.panel("Easing", default_closed=True):
            interpolation = tree.inputs.menu(
                "Interpolation",
                description="Shape of the easing curve, following Robert Penner's easing functions (easings.net)",
            )
            ease = tree.inputs.menu(
                "Ease",
                description="Apply the curve at the start (In), the end (Out) or both ends (In Out) of the transition",
            )
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
        value_1 = tree.outputs.float("Value")
        active = tree.outputs.boolean(
            "Active",
            description="Whether the input `Value` is between (and including) the lower and upper bounds",
        )
        end = tree.outputs.float("End")

        math_1 = start + length
        scene_time = g.SceneTime()
        menu_switch = g.MenuSwitch.float(
            time,
            {
                "Seconds": (
                    scene_time.o.seconds,
                    "Use the current time in `Seconds` from the active scene to aniamte the action",
                ),
                "Frames": (
                    scene_time.o.frame,
                    "Use the current time in `Frames` from the active scene",
                ),
                "Value": (value, "Use an arbitrary value to animate the action"),
            },
        )
        BetweenFloat(value=menu_switch.o.output, lower=start, upper=math_1) >> active
        map_range = menu_switch.o.output.map_range(start, math_1)
        switch = stagger.switch.float(
            map_range,
            StaggerValue(value=map_range, width=width, id=id, group_id=group_id),
        )
        (
            EaseValue(
                value=switch,
                interpolation=interpolation,
                ease=ease,
                overshoot=overshoot,
                period=period,
                amplitude=amplitude,
            )
            >> value_1
        )

        math_1 >> end

        time.default_value = "Seconds"
        interpolation.default_value = "Linear"
        ease.default_value = "In"


ASSET = AnimateAction

ASSET_METADATA = {
    "catalog_id": "85730213-4c2e-469f-b333-52ac53adf274",
}
