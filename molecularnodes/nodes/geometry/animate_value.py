# Node-group asset "Animate Value" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputFloat, InputInteger


class AnimateValue(AssetGeometryGroup):
    """
    Animate Value

    Parameters
    ----------
    smoother_step : InputBoolean
        Ease out and in from the min and max values
    clamped : InputBoolean
        Whether to clamp the interpolated value to the max
    frame_start : InputInteger
        Frame to start the animation on
    frame_end : InputInteger
        Frame to finish the animation on
    value_min : InputFloat
        Value to start animation from
    value_max : InputFloat
        Value to end animation at

    Inputs
    ------
    i.smoother_step : BooleanSocket
        Ease out and in from the min and max values
    i.clamped : BooleanSocket
        Whether to clamp the interpolated value to the max
    i.frame_start : IntegerSocket
        Frame to start the animation on
    i.frame_end : IntegerSocket
        Frame to finish the animation on
    i.value_min : FloatSocket
        Value to start animation from
    i.value_max : FloatSocket
        Value to end animation at

    Outputs
    -------
    o.value : FloatSocket
        Animated value that interpolates from min to max over frames
    """

    _name = "Animate Value"
    _asset_name = "Animate Value"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.animate_value"}

    class _Inputs(SocketAccessor):
        smoother_step: BooleanSocket
        """Ease out and in from the min and max values"""
        clamped: BooleanSocket
        """Whether to clamp the interpolated value to the max"""
        frame_start: IntegerSocket
        """Frame to start the animation on"""
        frame_end: IntegerSocket
        """Frame to finish the animation on"""
        value_min: FloatSocket
        """Value to start animation from"""
        value_max: FloatSocket
        """Value to end animation at"""

    class _Outputs(SocketAccessor):
        value: FloatSocket
        """Animated value that interpolates from min to max over frames"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        smoother_step: InputBoolean = False,
        clamped: InputBoolean = False,
        frame_start: InputInteger = 1,
        frame_end: InputInteger = 250,
        value_min: InputFloat = 0.0,
        value_max: InputFloat = 1.0,
    ):
        super().__init__(
            **{
                "Smoother Step": smoother_step,
                "Clamped": clamped,
                "Frame Start": frame_start,
                "Frame End": frame_end,
                "Value Min": value_min,
                "Value Max": value_max,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        smoother_step = tree.inputs.boolean(
            "Smoother Step",
            False,
            description="Ease out and in from the min and max values",
        )
        clamped = tree.inputs.boolean(
            "Clamped",
            False,
            description="Whether to clamp the interpolated value to the max",
        )
        with tree.inputs.panel("Frame"):
            frame_start = tree.inputs.integer(
                "Frame Start",
                1,
                description="Frame to start the animation on",
                min_value=1,
            )
            frame_end = tree.inputs.integer(
                "Frame End",
                250,
                description="Frame to finish the animation on",
                min_value=1,
            )
        with tree.inputs.panel("Value"):
            value_min = tree.inputs.float(
                "Value Min",
                0.0,
                description="Value to start animation from",
                min_value=-10_000.0,
                max_value=10_000.0,
            )
            value_max = tree.inputs.float(
                "Value Max",
                1.0,
                description="Value to end animation at",
                min_value=-10_000.0,
                max_value=10_000.0,
            )
        value = tree.outputs.float(
            "Value",
            description="Animated value that interpolates from min to max over frames",
        )

        scene_time = g.SceneTime()
        map_range = scene_time.o.frame.map_range(
            frame_start,
            frame_end,
            value_min,
            value_max,
            clamp=False,
            interpolation_type="SMOOTHERSTEP",
        )
        switch = clamped.switch.float(
            scene_time.o.frame.map_range(
                frame_start, frame_end, value_min, value_max, clamp=False
            ),
            scene_time.o.frame.map_range(frame_start, frame_end, value_min, value_max),
        )
        smoother_step.switch.float(switch, map_range) >> value


ASSET = AnimateValue

ASSET_METADATA = {
    "catalog_id": "85730213-4c2e-469f-b333-52ac53adf274",
}
