# Node-group asset "Animate Stagger" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
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
from nodebpy.types import InputBoolean, InputFloat, InputInteger, InputMenu
from .animate_ease import AnimateEase
from .chain_id import ChainID
from .ures_id import UResID


class AnimateStagger(AssetGeometryGroup):
    """
    Animate Stagger

    Parameters
    ----------
    frame : InputFloat
        Frame to evaluate the stagger at. Defaults to the current scene frame
    order : InputMenu | Literal["Residue", "Chain", "Atom", "Attribute"]
        What delays each point's start: its residue (`URes ID`), its chain, its atom index or the `Attribute` input
    attribute : InputFloat
        Per-point value multiplied by `Delay` when `Order` is `Attribute`, e.g. `res_id`, or a stored start frame with `Frame Start` 0 and `Delay` 1
    reverse : InputBoolean
        Stagger from the highest `Order` value to the lowest, so the last residue starts first
    frame_start : InputInteger
        Frame the first point starts on
    delay : InputFloat
        Frames between the starts of consecutive residues, chains, atoms or attribute values
    length : InputFloat
        Frames each point takes to go from 0 to 1. At 0 the factor snaps
    interpolation : InputMenu | Literal["Linear", "Sinusoidal", "Quadratic", "Cubic", "Quartic", "Quintic", "Exponential", "Circular", "Back", "Bounce", "Elastic"]
        Shape of each point's transition, see `Animate Ease`
    ease : InputMenu | Literal["In", "Out", "In Out"]
        Apply the curve at the start (In), the end (Out) or both ends (In Out) of each transition

    Inputs
    ------
    i.frame : FloatSocket
        Frame to evaluate the stagger at. Defaults to the current scene frame
    i.order : MenuSocket
        What delays each point's start: its residue (`URes ID`), its chain, its atom index or the `Attribute` input
    i.attribute : FloatSocket
        Per-point value multiplied by `Delay` when `Order` is `Attribute`, e.g. `res_id`, or a stored start frame with `Frame Start` 0 and `Delay` 1
    i.reverse : BooleanSocket
        Stagger from the highest `Order` value to the lowest, so the last residue starts first
    i.frame_start : IntegerSocket
        Frame the first point starts on
    i.delay : FloatSocket
        Frames between the starts of consecutive residues, chains, atoms or attribute values
    i.length : FloatSocket
        Frames each point takes to go from 0 to 1. At 0 the factor snaps
    i.interpolation : MenuSocket
        Shape of each point's transition, see `Animate Ease`
    i.ease : MenuSocket
        Apply the curve at the start (In), the end (Out) or both ends (In Out) of each transition

    Outputs
    -------
    o.factor : FloatSocket
        0 before each point's start, 1 once it has run for `Length` frames and eased in between. Feed it to `Set Color`, `Animate Reveal` or a symmetry node's `Factor`
    o.start : FloatSocket
        Frame each point starts its transition on
    """

    _name = "Animate Stagger"
    _asset_name = "Animate Stagger"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        frame: FloatSocket
        """Frame to evaluate the stagger at. Defaults to the current scene frame"""
        order: MenuSocket
        """What delays each point's start: its residue (`URes ID`), its chain, its atom index or the `Attribute` input"""
        attribute: FloatSocket
        """Per-point value multiplied by `Delay` when `Order` is `Attribute`, e.g. `res_id`, or a stored start frame with `Frame Start` 0 and `Delay` 1"""
        reverse: BooleanSocket
        """Stagger from the highest `Order` value to the lowest, so the last residue starts first"""
        frame_start: IntegerSocket
        """Frame the first point starts on"""
        delay: FloatSocket
        """Frames between the starts of consecutive residues, chains, atoms or attribute values"""
        length: FloatSocket
        """Frames each point takes to go from 0 to 1. At 0 the factor snaps"""
        interpolation: MenuSocket
        """Shape of each point's transition, see `Animate Ease`"""
        ease: MenuSocket
        """Apply the curve at the start (In), the end (Out) or both ends (In Out) of each transition"""

    class _Outputs(SocketAccessor):
        factor: FloatSocket
        """0 before each point's start, 1 once it has run for `Length` frames and eased in between. Feed it to `Set Color`, `Animate Reveal` or a symmetry node's `Factor`"""
        start: FloatSocket
        """Frame each point starts its transition on"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        frame: InputFloat = 0.0,
        order: InputMenu | Literal["Residue", "Chain", "Atom", "Attribute"] = "Residue",
        attribute: InputFloat = 0.0,
        reverse: InputBoolean = False,
        frame_start: InputInteger = 1,
        delay: InputFloat = 1.0,
        length: InputFloat = 25.0,
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
    ):
        super().__init__(
            **{
                "Frame": frame,
                "Order": order,
                "Attribute": attribute,
                "Reverse": reverse,
                "Frame Start": frame_start,
                "Delay": delay,
                "Length": length,
                "Interpolation": interpolation,
                "Ease": ease,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        frame = tree.inputs.float(
            "Frame",
            0.0,
            description="Frame to evaluate the stagger at. Defaults to the current scene frame",
            default_input="SCENE_FRAME",
        )
        order = tree.inputs.menu(
            "Order",
            description="What delays each point's start: its residue (`URes ID`), its chain, its atom index or the `Attribute` input",
        )
        attribute = tree.inputs.float(
            "Attribute",
            0.0,
            description="Per-point value multiplied by `Delay` when `Order` is `Attribute`, e.g. `res_id`, or a stored start frame with `Frame Start` 0 and `Delay` 1",
            hide_value=True,
        )
        reverse = tree.inputs.boolean(
            "Reverse",
            False,
            description="Stagger from the highest `Order` value to the lowest, so the last residue starts first",
        )
        with tree.inputs.panel("Timing"):
            frame_start = tree.inputs.integer(
                "Frame Start", 1, description="Frame the first point starts on"
            )
            delay = tree.inputs.float(
                "Delay",
                1.0,
                description="Frames between the starts of consecutive residues, chains, atoms or attribute values",
                min_value=0.0,
                max_value=10_000.0,
            )
            length = tree.inputs.float(
                "Length",
                25.0,
                description="Frames each point takes to go from 0 to 1. At 0 the factor snaps",
                min_value=0.0,
                max_value=10_000.0,
            )
        with tree.inputs.panel("Easing"):
            interpolation = tree.inputs.menu(
                "Interpolation",
                description="Shape of each point's transition, see `Animate Ease`",
            )
            ease = tree.inputs.menu(
                "Ease",
                description="Apply the curve at the start (In), the end (Out) or both ends (In Out) of each transition",
            )
        factor = tree.outputs.float(
            "Factor",
            description="0 before each point's start, 1 once it has run for `Length` frames and eased in between. Feed it to `Set Color`, `Animate Reveal` or a symmetry node's `Factor`",
        )
        start = tree.outputs.float(
            "Start", description="Frame each point starts its transition on"
        )

        with g.Frame("Per-point start frame"):
            _string = g.String(
                string="rank is the value each point is delayed by. Reverse flips it against the field maximum so the last rank starts first. start = Frame Start + rank * Delay."
            )
            menu_switch = g.MenuSwitch.float(
                order,
                {
                    "Residue": UResID().o.ures_id,
                    "Chain": ChainID(),
                    "Atom": g.Index(),
                    "Attribute": attribute,
                },
            )
            switch = reverse.switch.float(
                menu_switch.o.output,
                menu_switch.o.output.point.max() - menu_switch.o.output,
            )
            math_1 = frame_start + switch * delay
        with g.Frame("Progress"):
            _string_1 = g.String(
                string="Linear progress of each point over Length frames. Map Range returns 0 when its range is empty, so a Length of 0 falls back to a step at the start frame."
            )
            switch_1 = (length <= 0.0).switch.float(
                frame.map_range(math_1, math_1 + length, clamp=False),
                (frame >= math_1).switch.float(true=1.0),
            )
        AnimateEase(value=switch_1, interpolation=interpolation, ease=ease) >> factor

        math_1 >> start

        order.default_value = "Residue"
        interpolation.default_value = "Cubic"
        ease.default_value = "In Out"


ASSET = AnimateStagger

ASSET_METADATA = {
    "description": "Per-point 0..1 factor that starts later for each residue, chain, atom or attribute value",
    "catalog_id": "85730213-4c2e-469f-b333-52ac53adf274",
}
