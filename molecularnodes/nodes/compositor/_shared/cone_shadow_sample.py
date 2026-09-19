# Node group "Cone Shadow Sample" (CompositorNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import CompositorNodeTree
from nodebpy import TreeBuilder
from nodebpy import compositor as c
from nodebpy import geometry as g
from nodebpy.builder import CustomCompositorGroup, FloatSocket, SocketAccessor
from nodebpy.types import InputFloat


class ConeShadowSample(CustomCompositorGroup):
    """
    One sample of the conical shadow: whether the depth at a pixel offset is closer to the camera than the pixel by more than a threshold

    Parameters
    ----------
    depth : InputFloat
        Depth pass, in world units
    radius : InputFloat
        Distance of the sample from the pixel, in pixels
    angle : InputFloat
        Direction of the sample from the pixel
    threshold : InputFloat
        Depth gap, in world units, above which the sample casts shadow on the pixel

    Inputs
    ------
    i.depth : FloatSocket
        Depth pass, in world units
    i.radius : FloatSocket
        Distance of the sample from the pixel, in pixels
    i.angle : FloatSocket
        Direction of the sample from the pixel
    i.threshold : FloatSocket
        Depth gap, in world units, above which the sample casts shadow on the pixel

    Outputs
    -------
    o.occluded : FloatSocket
        Occluded
    """

    _name = "Cone Shadow Sample"
    _color_tag = "FILTER"
    _tree_properties = {
        "description": "One sample of the conical shadow: whether the depth at a pixel offset is closer to the camera than the pixel by more than a threshold"
    }

    class _Inputs(SocketAccessor):
        depth: FloatSocket
        """Depth pass, in world units"""
        radius: FloatSocket
        """Distance of the sample from the pixel, in pixels"""
        angle: FloatSocket
        """Direction of the sample from the pixel"""
        threshold: FloatSocket
        """Depth gap, in world units, above which the sample casts shadow on the pixel"""

    class _Outputs(SocketAccessor):
        occluded: FloatSocket
        """Occluded"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        depth: InputFloat = 0.0,
        radius: InputFloat = 0.0,
        angle: InputFloat = 0.0,
        threshold: InputFloat = 0.0,
    ):
        super().__init__(
            **{"Depth": depth, "Radius": radius, "Angle": angle, "Threshold": threshold}
        )

    def _build_group(self, tree: TreeBuilder[CompositorNodeTree]) -> None:
        depth = tree.inputs.float(
            "Depth", 0.0, description="Depth pass, in world units", hide_value=True
        )
        radius = tree.inputs.float(
            "Radius",
            0.0,
            description="Distance of the sample from the pixel, in pixels",
            min_value=0.0,
            max_value=10_000.0,
        )
        angle = tree.inputs.float(
            "Angle",
            0.0,
            description="Direction of the sample from the pixel",
            subtype="ANGLE",
        )
        threshold = tree.inputs.float(
            "Threshold",
            0.0,
            description="Depth gap, in world units, above which the sample casts shadow on the pixel",
            min_value=0.0,
            max_value=10_000.0,
        )
        occluded = tree.outputs.float("Occluded")

        math_1 = depth.min(10_000.0)
        translate = c.Translate(
            image=math_1,
            x=radius * angle.cos(),
            y=radius * angle.sin(),
            interpolation="Nearest",
            extension_x="Extend",
            extension_y="Extend",
        )
        math_2 = g.Math.greater_than(math_1 - translate, threshold)

        math_2 >> occluded
