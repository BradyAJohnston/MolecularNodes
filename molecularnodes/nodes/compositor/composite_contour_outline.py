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
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputFloat
from ._shared.depth_difference_sample import DepthDifferenceSample


class CompositeContourOutline(AssetCompositorGroup):
    """
    Contour outline opacity from the Depth pass after Goodsell's Illustrate: a 3x3 Laplacian of the depth with each neighbour's difference clamped to [Min Diff, Max Diff] in Angstrom, ramped from Low to High into an opacity. Lines sit on the nearer surface of every depth step and their weight follows the size of the step, unlike the binary Sobel lines of Composite Outline Mask

    Parameters
    ----------
    depth : InputFloat
        Depth pass from the Render Layers node, in world units
    low : InputFloat
        Laplacian value, in Angstrom, at which the line starts to appear
    high : InputFloat
        Laplacian value, in Angstrom, at which the line is fully opaque. A narrower range gives harder lines
    min_diff : InputFloat
        Smallest depth step between neighbouring pixels, in Angstrom, that contributes to a line
    max_diff : InputFloat
        Depth steps larger than this, in Angstrom, count as this much; a wider range emphasises larger features
    smooth : InputBoolean
        Average the opacity over the 3x3 neighbourhood where at least six of its nine pixels carry line signal, softening jagged lines
    world_scale : InputFloat
        World units per Angstrom, used to convert the Angstrom inputs. Molecular Nodes imports structures at 0.1 (1 nm per world unit)

    Inputs
    ------
    i.depth : FloatSocket
        Depth pass from the Render Layers node, in world units
    i.low : FloatSocket
        Laplacian value, in Angstrom, at which the line starts to appear
    i.high : FloatSocket
        Laplacian value, in Angstrom, at which the line is fully opaque. A narrower range gives harder lines
    i.min_diff : FloatSocket
        Smallest depth step between neighbouring pixels, in Angstrom, that contributes to a line
    i.max_diff : FloatSocket
        Depth steps larger than this, in Angstrom, count as this much; a wider range emphasises larger features
    i.smooth : BooleanSocket
        Average the opacity over the 3x3 neighbourhood where at least six of its nine pixels carry line signal, softening jagged lines
    i.world_scale : FloatSocket
        World units per Angstrom, used to convert the Angstrom inputs. Molecular Nodes imports structures at 0.1 (1 nm per world unit)

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
        low: FloatSocket
        """Laplacian value, in Angstrom, at which the line starts to appear"""
        high: FloatSocket
        """Laplacian value, in Angstrom, at which the line is fully opaque. A narrower range gives harder lines"""
        min_diff: FloatSocket
        """Smallest depth step between neighbouring pixels, in Angstrom, that contributes to a line"""
        max_diff: FloatSocket
        """Depth steps larger than this, in Angstrom, count as this much; a wider range emphasises larger features"""
        smooth: BooleanSocket
        """Average the opacity over the 3x3 neighbourhood where at least six of its nine pixels carry line signal, softening jagged lines"""
        world_scale: FloatSocket
        """World units per Angstrom, used to convert the Angstrom inputs. Molecular Nodes imports structures at 0.1 (1 nm per world unit)"""

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
        low: InputFloat = 3.0,
        high: InputFloat = 10.0,
        min_diff: InputFloat = 0.0,
        max_diff: InputFloat = 5.0,
        smooth: InputBoolean = True,
        world_scale: InputFloat = 0.1,
    ):
        super().__init__(
            **{
                "Depth": depth,
                "Low": low,
                "High": high,
                "Min Diff": min_diff,
                "Max Diff": max_diff,
                "Smooth": smooth,
                "World Scale": world_scale,
            }
        )

    def _build_group(self, tree: TreeBuilder[CompositorNodeTree]) -> None:
        depth = tree.inputs.float(
            "Depth",
            0.0,
            description="Depth pass from the Render Layers node, in world units",
            hide_value=True,
        )
        low = tree.inputs.float(
            "Low",
            3.0,
            description="Laplacian value, in Angstrom, at which the line starts to appear",
            min_value=0.0,
            max_value=10_000.0,
        )
        high = tree.inputs.float(
            "High",
            10.0,
            description="Laplacian value, in Angstrom, at which the line is fully opaque. A narrower range gives harder lines",
            min_value=0.0,
            max_value=10_000.0,
        )
        min_diff = tree.inputs.float(
            "Min Diff",
            0.0,
            description="Smallest depth step between neighbouring pixels, in Angstrom, that contributes to a line",
            min_value=0.0,
            max_value=10_000.0,
        )
        max_diff = tree.inputs.float(
            "Max Diff",
            5.0,
            description="Depth steps larger than this, in Angstrom, count as this much; a wider range emphasises larger features",
            min_value=0.0,
            max_value=10_000.0,
        )
        smooth = tree.inputs.boolean(
            "Smooth",
            True,
            description="Average the opacity over the 3x3 neighbourhood where at least six of its nine pixels carry line signal, softening jagged lines",
        )
        world_scale = tree.inputs.float(
            "World Scale",
            0.1,
            description="World units per Angstrom, used to convert the Angstrom inputs. Molecular Nodes imports structures at 0.1 (1 nm per world unit)",
            min_value=0.0,
            max_value=10_000.0,
        )
        opacity = tree.outputs.float("Opacity")

        with c.Frame("Laplacian"):
            math_1 = min_diff * world_scale
            math_2 = max_diff * world_scale
            math_3 = DepthDifferenceSample(
                depth=depth, x=-1.0, y=-1.0, min=math_1, max=math_2, weight=0.8
            ).o.value + DepthDifferenceSample(
                depth=depth, y=-1.0, min=math_1, max=math_2
            )
            math_4 = math_3 + DepthDifferenceSample(
                depth=depth, x=1.0, y=-1.0, min=math_1, max=math_2, weight=0.8
            )
            math_5 = (
                math_4
                + DepthDifferenceSample(depth=depth, x=-1.0, min=math_1, max=math_2)
                + DepthDifferenceSample(depth=depth, x=1.0, min=math_1, max=math_2)
            )
            math_6 = math_5 + DepthDifferenceSample(
                depth=depth, x=-1.0, y=1.0, min=math_1, max=math_2, weight=0.8
            )
            math_7 = (
                math_6
                + DepthDifferenceSample(depth=depth, y=1.0, min=math_1, max=math_2)
                + DepthDifferenceSample(
                    depth=depth, x=1.0, y=1.0, min=math_1, max=math_2, weight=0.8
                )
            )
            map_range = (math_7 / world_scale).map_range(low, high)
            _string = g.String(
                string="Each of the eight neighbours contributes how much farther it is than the centre pixel, clamped to [Min Diff, Max Diff], weighted 1.0 on the edges and 0.8 on the corners. Nearer neighbours contribute nothing, so the line lands on the nearer side of a depth step."
            )
        with c.Frame("Smoothing"):
            math_8 = g.Math.greater_than(
                c.Blur(
                    image=g.Math.greater_than(map_range, 0.0),
                    size=(1.0, 1.0),
                    type="Flat",
                ),
                0.6,
            )
            mix = math_8.o.value.mix.float(
                map_range, c.Blur(image=map_range, size=(1.0, 1.0), type="Flat")
            )
            mix_1 = g.Mix.float(smooth, map_range, mix)
            _string_1 = g.String(
                string="Illustrate averages the 3x3 neighbourhood when at least six of its pixels carry signal. A flat blur of radius one is that average, and the same blur of the signal mask is the fraction of pixels that carry it."
            )

        mix_1 >> opacity


ASSET = CompositeContourOutline

ASSET_METADATA = {
    "description": "Depth Laplacian contour outlines after Goodsell's Illustrate",
    "catalog_id": "441e6ca5-e514-4e77-a3cd-25fc1a2e08ae",
}
