# Node group "Depth Difference Sample" (CompositorNodeTree), dumped by nodebpy.assets.dump_library.
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


class DepthDifferenceSample(CustomCompositorGroup):
    """
    One term of the depth Laplacian: how much farther the pixel at an offset is than the centre, clamped to [Min, Max] and weighted

    Parameters
    ----------
    depth : InputFloat
        Depth pass, in world units
    x : InputFloat
        Offset of the sample in pixels
    y : InputFloat
        Offset of the sample in pixels
    min : InputFloat
        Smallest depth difference counted, in world units
    max : InputFloat
        Largest depth difference counted, in world units
    weight : InputFloat
        Kernel weight

    Inputs
    ------
    i.depth : FloatSocket
        Depth pass, in world units
    i.x : FloatSocket
        Offset of the sample in pixels
    i.y : FloatSocket
        Offset of the sample in pixels
    i.min : FloatSocket
        Smallest depth difference counted, in world units
    i.max : FloatSocket
        Largest depth difference counted, in world units
    i.weight : FloatSocket
        Kernel weight

    Outputs
    -------
    o.value : FloatSocket
        Value
    """

    _name = "Depth Difference Sample"
    _color_tag = "FILTER"
    _tree_properties = {
        "description": "One term of the depth Laplacian: how much farther the pixel at an offset is than the centre, clamped to [Min, Max] and weighted"
    }

    class _Inputs(SocketAccessor):
        depth: FloatSocket
        """Depth pass, in world units"""
        x: FloatSocket
        """Offset of the sample in pixels"""
        y: FloatSocket
        """Offset of the sample in pixels"""
        min: FloatSocket
        """Smallest depth difference counted, in world units"""
        max: FloatSocket
        """Largest depth difference counted, in world units"""
        weight: FloatSocket
        """Kernel weight"""

    class _Outputs(SocketAccessor):
        value: FloatSocket
        """Value"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        depth: InputFloat = 0.0,
        x: InputFloat = 0.0,
        y: InputFloat = 0.0,
        min: InputFloat = 0.0,
        max: InputFloat = 1.0,
        weight: InputFloat = 1.0,
    ):
        super().__init__(
            **{"Depth": depth, "X": x, "Y": y, "Min": min, "Max": max, "Weight": weight}
        )

    def _build_group(self, tree: TreeBuilder[CompositorNodeTree]) -> None:
        depth = tree.inputs.float(
            "Depth", 0.0, description="Depth pass, in world units", hide_value=True
        )
        x = tree.inputs.float("X", 0.0, description="Offset of the sample in pixels")
        y = tree.inputs.float("Y", 0.0, description="Offset of the sample in pixels")
        min = tree.inputs.float(
            "Min", 0.0, description="Smallest depth difference counted, in world units"
        )
        max = tree.inputs.float(
            "Max", 1.0, description="Largest depth difference counted, in world units"
        )
        weight = tree.inputs.float("Weight", 1.0, description="Kernel weight")
        value = tree.outputs.float("Value")

        math_1 = depth.min(10_000.0)
        translate = c.Translate(
            image=math_1,
            x=x,
            y=y,
            interpolation="Nearest",
            extension_x="Extend",
            extension_y="Extend",
        )
        g.Math.subtract(translate, math_1).o.value.clamp(min, max) * weight >> value
