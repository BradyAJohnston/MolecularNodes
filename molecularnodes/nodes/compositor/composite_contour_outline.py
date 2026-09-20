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
from ._shared.angstromtoworld import AngstromToWorld2


class CompositeContourOutline(AssetCompositorGroup):
    """
    Contour outline opacity from the Depth pass after Goodsell's Illustrate: a 3x3 Laplacian of the depth with each neighbour's difference clamped to [Min Diff, Max Diff] in Angstrom, ramped from Low to High into an opacity. Lines sit on the nearer surface of every depth step and their weight follows the size of the step, unlike the binary Sobel lines of Composite Outline Mask

    Parameters
    ----------
    depth : InputFloat
        Depth pass from the Render Layers node, in world units
    high : InputFloat
        Laplacian value, in Angstrom, at which the line is fully opaque. A narrower range gives harder lines
    smooth : InputBoolean
        Average the opacity over the 3x3 neighbourhood where at least six of its nine pixels carry line signal, softening jagged lines

    Inputs
    ------
    i.depth : FloatSocket
        Depth pass from the Render Layers node, in world units
    i.high : FloatSocket
        Laplacian value, in Angstrom, at which the line is fully opaque. A narrower range gives harder lines
    i.smooth : BooleanSocket
        Average the opacity over the 3x3 neighbourhood where at least six of its nine pixels carry line signal, softening jagged lines

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
        high: FloatSocket
        """Laplacian value, in Angstrom, at which the line is fully opaque. A narrower range gives harder lines"""
        smooth: BooleanSocket
        """Average the opacity over the 3x3 neighbourhood where at least six of its nine pixels carry line signal, softening jagged lines"""

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
        high: InputFloat = 10.0,
        smooth: InputBoolean = True,
    ):
        super().__init__(**{"Depth": depth, "High": high, "Smooth": smooth})

    def _build_group(self, tree: TreeBuilder[CompositorNodeTree]) -> None:
        depth = tree.inputs.float(
            "Depth",
            0.0,
            description="Depth pass from the Render Layers node, in world units",
            hide_value=True,
        )
        high = tree.inputs.float(
            "High",
            10.0,
            description="Laplacian value, in Angstrom, at which the line is fully opaque. A narrower range gives harder lines",
            min_value=0.0,
            max_value=10_000.0,
        )
        smooth = tree.inputs.boolean(
            "Smooth",
            True,
            description="Average the opacity over the 3x3 neighbourhood where at least six of its nine pixels carry line signal, softening jagged lines",
        )
        opacity = tree.outputs.float("Opacity")

        with c.Frame("Laplacian"):
            map_range = AngstromToWorld2(
                angstrom=c.Filter(image=depth, type="Laplace")
            ).o.world.map_range(from_max=high)
        with c.Frame("Smoothing"):
            math_1 = g.Math.greater_than(
                c.Blur(
                    image=g.Math.greater_than(map_range, 0.0),
                    size=(1.0, 1.0),
                    type="Flat",
                ),
                0.6,
            )
            mix = math_1.o.value.mix.float(
                map_range, c.Blur(image=map_range, size=(1.0, 1.0), type="Flat")
            )
            mix_1 = g.Mix.float(smooth, map_range, mix)
            _string = g.String(
                string="Illustrate averages the 3x3 neighbourhood when at least six of its pixels carry signal. A flat blur of radius one is that average, and the same blur of the signal mask is the fraction of pixels that carry it."
            )

        mix_1 >> opacity


ASSET = CompositeContourOutline

ASSET_METADATA = {
    "description": "Depth Laplacian contour outlines after Goodsell's Illustrate",
    "catalog_id": "441e6ca5-e514-4e77-a3cd-25fc1a2e08ae",
}
