# Node-group asset "Depth Outline" (CompositorNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import CompositorNodeTree
from nodebpy import TreeBuilder
from nodebpy import compositor as c
from nodebpy.builder import (
    AssetCompositorGroup,
    BooleanSocket,
    FloatSocket,
    IntegerSocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputBoolean, InputFloat, InputInteger, InputMenu, InputVector
from ._shared.depth_difference_sample import DepthDifferenceSample


class DepthOutline(AssetCompositorGroup):
    """
    Line mask from the Depth pass, by one of two methods: Sobel, a thresholded Sobel edge detector on depth, or Goodsell, Illustrate's 3x3 Laplacian with each neighbour's step clamped to [Min Diff, Max Diff] and ramped from Low to High into a soft opacity. Optionally adds lines where the normal changes sharply; the lines are anti-aliased and thickened to Size pixels

    Parameters
    ----------
    depth : InputFloat
        Depth pass from the Render Layers node, in world units
    method : InputMenu | Literal["Sobel", "Goodsell"]
        Sobel draws hard lines wherever the depth gradient passes Threshold; Goodsell draws soft lines on the nearer side of each depth step, weighted by the size of the step, as in Illustrate
    size : InputInteger
        Line thickness in pixels
    threshold : InputFloat
        Depth gradient, in world units like the depth pass, above which a line is drawn, for the Sobel method
    low : InputFloat
        Weighted sum of the clamped depth steps, in world units, at which a line starts to appear
    high : InputFloat
        Weighted sum of the clamped depth steps, in world units, at which a line is fully opaque. A narrower range from Low gives harder lines
    min_diff : InputFloat
        Smallest depth step between neighbouring pixels, in world units, that contributes to a line
    max_diff : InputFloat
        Depth steps larger than this, in world units, count as this much; a wider range emphasises larger features
    smooth : InputBoolean
        Average the opacity over the 3x3 neighbourhood where at least six of its nine pixels carry line signal, softening jagged lines
    normals : InputBoolean
        Also draw lines where the surface normal changes sharply, catching creases between surfaces at the same depth
    normals_pass : InputVector
        Normal pass from the Render Layers node
    normals_threshold : InputFloat
        Change in normal between neighbouring pixels above which a line is drawn, when Normals is enabled

    Inputs
    ------
    i.depth : FloatSocket
        Depth pass from the Render Layers node, in world units
    i.method : MenuSocket
        Sobel draws hard lines wherever the depth gradient passes Threshold; Goodsell draws soft lines on the nearer side of each depth step, weighted by the size of the step, as in Illustrate
    i.size : IntegerSocket
        Line thickness in pixels
    i.threshold : FloatSocket
        Depth gradient, in world units like the depth pass, above which a line is drawn, for the Sobel method
    i.low : FloatSocket
        Weighted sum of the clamped depth steps, in world units, at which a line starts to appear
    i.high : FloatSocket
        Weighted sum of the clamped depth steps, in world units, at which a line is fully opaque. A narrower range from Low gives harder lines
    i.min_diff : FloatSocket
        Smallest depth step between neighbouring pixels, in world units, that contributes to a line
    i.max_diff : FloatSocket
        Depth steps larger than this, in world units, count as this much; a wider range emphasises larger features
    i.smooth : BooleanSocket
        Average the opacity over the 3x3 neighbourhood where at least six of its nine pixels carry line signal, softening jagged lines
    i.normals : BooleanSocket
        Also draw lines where the surface normal changes sharply, catching creases between surfaces at the same depth
    i.normals_pass : VectorSocket
        Normal pass from the Render Layers node
    i.normals_threshold : FloatSocket
        Change in normal between neighbouring pixels above which a line is drawn, when Normals is enabled

    Outputs
    -------
    o.mask : FloatSocket
        Mask
    """

    _name = "Depth Outline"
    _asset_name = "Depth Outline"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "FILTER"
    _tree_properties = {
        "description": "Line mask from the Depth pass, by one of two methods: Sobel, a thresholded Sobel edge detector on depth, or Goodsell, Illustrate's 3x3 Laplacian with each neighbour's step clamped to [Min Diff, Max Diff] and ramped from Low to High into a soft opacity. Optionally adds lines where the normal changes sharply; the lines are anti-aliased and thickened to Size pixels"
    }

    class _Inputs(SocketAccessor):
        depth: FloatSocket
        """Depth pass from the Render Layers node, in world units"""
        method: MenuSocket
        """Sobel draws hard lines wherever the depth gradient passes Threshold; Goodsell draws soft lines on the nearer side of each depth step, weighted by the size of the step, as in Illustrate"""
        size: IntegerSocket
        """Line thickness in pixels"""
        threshold: FloatSocket
        """Depth gradient, in world units like the depth pass, above which a line is drawn, for the Sobel method"""
        low: FloatSocket
        """Weighted sum of the clamped depth steps, in world units, at which a line starts to appear"""
        high: FloatSocket
        """Weighted sum of the clamped depth steps, in world units, at which a line is fully opaque. A narrower range from Low gives harder lines"""
        min_diff: FloatSocket
        """Smallest depth step between neighbouring pixels, in world units, that contributes to a line"""
        max_diff: FloatSocket
        """Depth steps larger than this, in world units, count as this much; a wider range emphasises larger features"""
        smooth: BooleanSocket
        """Average the opacity over the 3x3 neighbourhood where at least six of its nine pixels carry line signal, softening jagged lines"""
        normals: BooleanSocket
        """Also draw lines where the surface normal changes sharply, catching creases between surfaces at the same depth"""
        normals_pass: VectorSocket
        """Normal pass from the Render Layers node"""
        normals_threshold: FloatSocket
        """Change in normal between neighbouring pixels above which a line is drawn, when Normals is enabled"""

    class _Outputs(SocketAccessor):
        mask: FloatSocket
        """Mask"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        depth: InputFloat = 0.0,
        method: InputMenu | Literal["Sobel", "Goodsell"] = "Sobel",
        size: InputInteger = 1,
        threshold: InputFloat = 6.0,
        low: InputFloat = 0.3,
        high: InputFloat = 1.0,
        min_diff: InputFloat = 0.0,
        max_diff: InputFloat = 0.5,
        smooth: InputBoolean = True,
        normals: InputBoolean = False,
        normals_pass: InputVector = None,
        normals_threshold: InputFloat = 3.0,
    ):
        super().__init__(
            **{
                "Depth": depth,
                "Method": method,
                "Size": size,
                "Threshold": threshold,
                "Low": low,
                "High": high,
                "Min Diff": min_diff,
                "Max Diff": max_diff,
                "Smooth": smooth,
                "Normals": normals,
                "Normals Pass": normals_pass,
                "Normals Threshold": normals_threshold,
            }
        )

    def _build_group(self, tree: TreeBuilder[CompositorNodeTree]) -> None:
        depth = tree.inputs.float(
            "Depth",
            0.0,
            description="Depth pass from the Render Layers node, in world units",
            hide_value=True,
        )
        method = tree.inputs.menu(
            "Method",
            "Sobel",
            description="Sobel draws hard lines wherever the depth gradient passes Threshold; Goodsell draws soft lines on the nearer side of each depth step, weighted by the size of the step, as in Illustrate",
        )
        size = tree.inputs.integer(
            "Size", 1, description="Line thickness in pixels", min_value=1, max_value=20
        )
        threshold = tree.inputs.float(
            "Threshold",
            6.0,
            description="Depth gradient, in world units like the depth pass, above which a line is drawn, for the Sobel method",
            min_value=0.0,
            max_value=10_000.0,
        )
        with tree.inputs.panel("Goodsell", default_closed=True):
            low = tree.inputs.float(
                "Low",
                0.3,
                description="Weighted sum of the clamped depth steps, in world units, at which a line starts to appear",
                min_value=0.0,
                max_value=10_000.0,
            )
            high = tree.inputs.float(
                "High",
                1.0,
                description="Weighted sum of the clamped depth steps, in world units, at which a line is fully opaque. A narrower range from Low gives harder lines",
                min_value=0.0,
                max_value=10_000.0,
            )
            min_diff = tree.inputs.float(
                "Min Diff",
                0.0,
                description="Smallest depth step between neighbouring pixels, in world units, that contributes to a line",
                min_value=0.0,
                max_value=10_000.0,
            )
            max_diff = tree.inputs.float(
                "Max Diff",
                0.5,
                description="Depth steps larger than this, in world units, count as this much; a wider range emphasises larger features",
                min_value=0.0,
                max_value=10_000.0,
            )
            smooth = tree.inputs.boolean(
                "Smooth",
                True,
                description="Average the opacity over the 3x3 neighbourhood where at least six of its nine pixels carry line signal, softening jagged lines",
            )
        with tree.inputs.panel("Normals", default_closed=True):
            normals = tree.inputs.boolean(
                "Normals",
                False,
                description="Also draw lines where the surface normal changes sharply, catching creases between surfaces at the same depth",
                is_panel_toggle=True,
            )
            normals_pass = tree.inputs.vector(
                "Normals Pass",
                (0.0, 0.0, 0.0),
                description="Normal pass from the Render Layers node",
                hide_value=True,
            )
            normals_threshold = tree.inputs.float(
                "Normals Threshold",
                3.0,
                description="Change in normal between neighbouring pixels above which a line is drawn, when Normals is enabled",
                min_value=0.0,
                max_value=10_000.0,
            )
        mask = tree.outputs.float("Mask")

        with c.Frame("Goodsell"):
            depth_difference_sample = DepthDifferenceSample(
                depth=depth, x=-1.0, y=-1.0, min=min_diff, max=max_diff, weight=0.8
            )
            depth_difference_sample_1 = DepthDifferenceSample(
                depth=depth, x=1.0, y=-1.0, min=min_diff, max=max_diff, weight=0.8
            )
            depth_difference_sample_2 = DepthDifferenceSample(
                depth=depth, x=-1.0, y=1.0, min=min_diff, max=max_diff, weight=0.8
            )
            math_1 = depth_difference_sample.o.value + DepthDifferenceSample(
                depth=depth, y=-1.0, min=min_diff, max=max_diff
            )
            math_2 = (
                math_1
                + depth_difference_sample_1
                + DepthDifferenceSample(depth=depth, x=-1.0, min=min_diff, max=max_diff)
            )
            math_3 = (
                math_2
                + DepthDifferenceSample(depth=depth, x=1.0, min=min_diff, max=max_diff)
                + depth_difference_sample_2
            )
            math_4 = (
                math_3
                + DepthDifferenceSample(depth=depth, y=1.0, min=min_diff, max=max_diff)
                + DepthDifferenceSample(
                    depth=depth, x=1.0, y=1.0, min=min_diff, max=max_diff, weight=0.8
                )
            )
            map_range = math_4.map_range(low, high)
            math_5 = c.Math.greater_than(
                c.Blur(
                    image=c.Math.greater_than(map_range, 0.0),
                    size=(1.0, 1.0),
                    type="Flat",
                ),
                0.6,
            )
            mix = math_5.o.value.mix.float(
                map_range, c.Blur(image=map_range, size=(1.0, 1.0), type="Flat")
            )
            mix_1 = c.Mix.float(smooth, map_range, mix)
            _string = c.String(
                "Each of the eight neighbours contributes how much farther it is than the centre pixel, clamped to [Min Diff, Max Diff], weighted 1.0 on the edges and 0.8 on the corners, so the line lands on the nearer side of a depth step. Illustrate then averages the 3x3 neighbourhood when at least six of its pixels carry signal: a flat blur of radius one is that average, and the same blur of the signal mask is the fraction of pixels that carry it."
            )
        with c.Frame("Normal edges"):
            math_6 = c.Math.greater_than(
                c.VectorMath.length(c.Filter(image=normals_pass, type="Sobel")).o.value,
                normals_threshold,
            )
            math_7 = math_6.o.value * normals
        with c.Frame("Sobel"):
            math_8 = c.Math.greater_than(
                c.Filter(image=depth.min(10_000.0), type="Sobel"), threshold
            )
            _string_1 = c.String(
                "The depth pass is unbounded where nothing was rendered, so it is clipped before the Sobel filter to keep the gradient finite. Silhouettes against the background still exceed any sensible threshold, so they are outlined."
            )
        with c.Frame("Line weight"):
            math_9 = c.MenuSwitch.float(
                method, {"Sobel": math_8, "Goodsell": mix_1.o.result_float}
            ).o.output.max(math_7)
            anti_aliasing = c.AntiAliasing(
                c.DilateErode(
                    mask=c.AntiAliasing(math_9),
                    size=(size - 1.0).max(0.0),
                    type="Distance",
                )
            )

        anti_aliasing >> mask


ASSET = DepthOutline

ASSET_METADATA = {
    "description": "Line mask from the Depth and Normal render passes",
    "catalog_id": "874a72a6-7cf2-456a-ab07-e79f5a1e4294",
}
