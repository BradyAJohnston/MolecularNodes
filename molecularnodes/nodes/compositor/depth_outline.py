# Node-group asset "Depth Outline" (CompositorNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import CompositorNodeTree
from nodebpy import TreeBuilder
from nodebpy import compositor as c
from nodebpy.builder import (
    AssetCompositorGroup,
    BooleanSocket,
    FloatSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputBoolean, InputFloat, InputInteger, InputVector


class DepthOutline(AssetCompositorGroup):
    """
    Line mask from the Depth and Normal render passes: a Sobel edge detector on depth, optionally combined with one on normals, anti-aliased and dilated to the requested pixel width

    Parameters
    ----------
    depth : InputFloat
        Depth pass from the Render Layers node, in world units
    size : InputInteger
        Line thickness in pixels
    depth_threshold : InputFloat
        Depth difference between neighbouring pixels, in Angstrom, above which a line is drawn
    socket_3 : InputBoolean
        Also draw lines where the surface normal changes sharply, catching creases between surfaces at the same depth
    socket_1 : InputVector
        Normal pass from the Render Layers node
    normals_threshold : InputFloat
        Change in normal between neighbouring pixels above which a line is drawn, when Use Normals is enabled

    Inputs
    ------
    i.depth : FloatSocket
        Depth pass from the Render Layers node, in world units
    i.size : IntegerSocket
        Line thickness in pixels
    i.depth_threshold : FloatSocket
        Depth difference between neighbouring pixels, in Angstrom, above which a line is drawn
    i.socket_3 : BooleanSocket
        Also draw lines where the surface normal changes sharply, catching creases between surfaces at the same depth
    i.socket_1 : VectorSocket
        Normal pass from the Render Layers node
    i.normals_threshold : FloatSocket
        Change in normal between neighbouring pixels above which a line is drawn, when Use Normals is enabled

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
        "description": "Line mask from the Depth and Normal render passes: a Sobel edge detector on depth, optionally combined with one on normals, anti-aliased and dilated to the requested pixel width"
    }

    class _Inputs(SocketAccessor):
        depth: FloatSocket
        """Depth pass from the Render Layers node, in world units"""
        size: IntegerSocket
        """Line thickness in pixels"""
        depth_threshold: FloatSocket
        """Depth difference between neighbouring pixels, in Angstrom, above which a line is drawn"""
        socket_3: BooleanSocket
        """Also draw lines where the surface normal changes sharply, catching creases between surfaces at the same depth"""
        socket_1: VectorSocket
        """Normal pass from the Render Layers node"""
        normals_threshold: FloatSocket
        """Change in normal between neighbouring pixels above which a line is drawn, when Use Normals is enabled"""

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
        size: InputInteger = 1,
        depth_threshold: InputFloat = 6.0,
        socket_3: InputBoolean = False,
        socket_1: InputVector = None,
        normals_threshold: InputFloat = 3.0,
    ):
        super().__init__(
            **{
                "Depth": depth,
                "Size": size,
                "Depth Threshold": depth_threshold,
                "Normals Threshold": normals_threshold,
            },
            _named_links=[("Normals", socket_3), ("Normals", socket_1)],
        )

    def _build_group(self, tree: TreeBuilder[CompositorNodeTree]) -> None:
        depth = tree.inputs.float(
            "Depth",
            0.0,
            description="Depth pass from the Render Layers node, in world units",
            hide_value=True,
        )
        size = tree.inputs.integer(
            "Size", 1, description="Line thickness in pixels", min_value=1, max_value=20
        )
        depth_threshold = tree.inputs.float(
            "Depth Threshold",
            6.0,
            description="Depth difference between neighbouring pixels, in Angstrom, above which a line is drawn",
            min_value=0.0,
            max_value=10_000.0,
        )
        with tree.inputs.panel("Normals", default_closed=True):
            normals = tree.inputs.boolean(
                "Normals",
                False,
                description="Also draw lines where the surface normal changes sharply, catching creases between surfaces at the same depth",
                is_panel_toggle=True,
            )
            normals_1 = tree.inputs.vector(
                "Normals",
                (0.0, 0.0, 0.0),
                description="Normal pass from the Render Layers node",
                hide_value=True,
            )
            normals_threshold = tree.inputs.float(
                "Normals Threshold",
                3.0,
                description="Change in normal between neighbouring pixels above which a line is drawn, when Use Normals is enabled",
                min_value=0.0,
                max_value=10_000.0,
            )
        mask = tree.outputs.float("Mask")

        with c.Frame("Depth edges"):
            math_1 = c.Math.greater_than(
                c.Filter(image=depth.min(10_000.0), type="Sobel"), depth_threshold
            )
            _string = c.String(
                "The depth pass is unbounded where nothing was rendered, so it is clipped before the Sobel filter to keep the gradient finite. Silhouettes against the background still exceed any sensible threshold, so they are outlined."
            )
        with c.Frame("Normal edges"):
            math_2 = c.Math.greater_than(
                c.VectorMath.length(c.Filter(image=normals_1, type="Sobel")).o.value,
                normals_threshold,
            )
            math_3 = math_1.o.value.max(math_2.o.value * normals)
        with c.Frame("Line weight"):
            anti_aliasing = c.AntiAliasing(
                c.DilateErode(
                    mask=c.AntiAliasing(math_3),
                    size=(size - 1.0).max(0.0),
                    type="Distance",
                )
            )
        _viewer = c.Viewer(math_1)

        anti_aliasing >> mask


ASSET = DepthOutline

ASSET_METADATA = {
    "description": "Line mask from the Depth and Normal render passes",
    "catalog_id": "874a72a6-7cf2-456a-ab07-e79f5a1e4294",
}
