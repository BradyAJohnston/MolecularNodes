# Node-group asset "Composite Contour Outline" (CompositorNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import CompositorNodeTree
from nodebpy import TreeBuilder
from nodebpy import compositor as c
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetCompositorGroup,
    BooleanSocket,
    FloatSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputFloat, InputInteger


class CompositeContourOutline(AssetCompositorGroup):
    """
    Contour outline opacity from the Depth pass after Goodsell's Illustrate: a 3x3 Laplacian of the depth with each neighbour's difference clamped to [Min Diff, Max Diff] in Angstrom, ramped from Low to High into an opacity. Lines sit on the nearer surface of every depth step and their weight follows the size of the step, unlike the binary Sobel lines of Composite Outline Mask

    Parameters
    ----------
    depth : InputFloat
        Depth pass from the Render Layers node, in world units
    smooth : InputBoolean
        Average the opacity over the 3x3 neighbourhood where at least six of its nine pixels carry line signal, softening jagged lines
    value : InputFloat
        Value
    size : InputInteger
        The size of dilation/erosion in pixels. Positive values dilates and negative values erodes

    Inputs
    ------
    i.depth : FloatSocket
        Depth pass from the Render Layers node, in world units
    i.smooth : BooleanSocket
        Average the opacity over the 3x3 neighbourhood where at least six of its nine pixels carry line signal, softening jagged lines
    i.value : FloatSocket
        Value
    i.size : IntegerSocket
        The size of dilation/erosion in pixels. Positive values dilates and negative values erodes

    Outputs
    -------
    o.opacity : FloatSocket
        Opacity
    """

    _name = "Composite Contour Outline"
    _asset_name = "Composite Contour Outline"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "FILTER"
    _tree_properties = {
        "description": "Contour outline opacity from the Depth pass after Goodsell's Illustrate: a 3x3 Laplacian of the depth with each neighbour's difference clamped to [Min Diff, Max Diff] in Angstrom, ramped from Low to High into an opacity. Lines sit on the nearer surface of every depth step and their weight follows the size of the step, unlike the binary Sobel lines of Composite Outline Mask"
    }

    class _Inputs(SocketAccessor):
        depth: FloatSocket
        """Depth pass from the Render Layers node, in world units"""
        smooth: BooleanSocket
        """Average the opacity over the 3x3 neighbourhood where at least six of its nine pixels carry line signal, softening jagged lines"""
        value: FloatSocket
        """Value"""
        size: IntegerSocket
        """The size of dilation/erosion in pixels. Positive values dilates and negative values erodes"""

    class _Outputs(SocketAccessor):
        opacity: FloatSocket
        """Opacity"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        depth: InputFloat = 0.0,
        smooth: InputBoolean = True,
        value: InputFloat = 4.59,
        size: InputInteger = 0,
    ):
        super().__init__(
            **{"Depth": depth, "Smooth": smooth, "Value": value, "Size": size}
        )

    def _build_group(self, tree: TreeBuilder[CompositorNodeTree]) -> None:
        depth = tree.inputs.float(
            "Depth",
            0.0,
            description="Depth pass from the Render Layers node, in world units",
            hide_value=True,
        )
        smooth = tree.inputs.boolean(
            "Smooth",
            True,
            description="Average the opacity over the 3x3 neighbourhood where at least six of its nine pixels carry line signal, softening jagged lines",
        )
        value = tree.inputs.float(
            "Value", 4.59, min_value=-10_000.0, max_value=10_000.0
        )
        size = tree.inputs.integer(
            "Size",
            0,
            description="The size of dilation/erosion in pixels. Positive values dilates and negative values erodes",
            subtype="PIXEL",
        )
        opacity = tree.outputs.float("Opacity")

        with c.Frame("Laplacian"):
            math_1 = g.Math.greater_than(c.Filter(image=depth, type="Sobel"), value)
        with c.Frame("Smoothing"):
            math_2 = g.Math.greater_than(
                c.Blur(
                    image=g.Math.greater_than(math_1, 0.0), size=(1.0, 1.0), type="Flat"
                ),
                0.6,
            )
            mix = g.Mix.float(
                smooth,
                math_1,
                math_2.o.value.mix.float(
                    math_1, c.Blur(image=math_1, size=(1.0, 1.0), type="Flat")
                ),
            )
            dilate_erode = c.DilateErode(
                mask=c.AntiAliasing(image=mix.o.result_float, threshold=0.2),
                size=size,
                type="Distance",
            )
            _string = g.String(
                string="Illustrate averages the 3x3 neighbourhood when at least six of its pixels carry signal. A flat blur of radius one is that average, and the same blur of the signal mask is the fraction of pixels that carry it."
            )

        dilate_erode >> opacity


ASSET = CompositeContourOutline

ASSET_METADATA = {
    "description": "Depth Laplacian contour outlines after Goodsell's Illustrate",
    "catalog_id": "441e6ca5-e514-4e77-a3cd-25fc1a2e08ae",
}
