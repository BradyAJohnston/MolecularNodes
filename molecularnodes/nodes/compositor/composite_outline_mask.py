# Node-group asset "Composite Outline Mask" (CompositorNodeTree), dumped by nodebpy.assets.dump_library.
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
    VectorSocket,
)
from nodebpy.types import InputBoolean, InputFloat, InputInteger, InputVector


class CompositeOutlineMask(AssetCompositorGroup):
    """
    Line mask from the Depth and Normal render passes: a Sobel edge detector on depth, optionally combined with one on normals, anti-aliased and dilated to the requested pixel width

    Parameters
    ----------
    depth : InputFloat
        Depth pass from the Render Layers node, in world units
    normal : InputVector
        Normal pass from the Render Layers node
    depth_threshold : InputFloat
        Depth difference between neighbouring pixels, in Angstrom, above which a line is drawn
    use_normals : InputBoolean
        Also draw lines where the surface normal changes sharply, catching creases between surfaces at the same depth
    normal_threshold : InputFloat
        Change in normal between neighbouring pixels above which a line is drawn, when Use Normals is enabled
    size : InputInteger
        Line thickness in pixels
    world_scale : InputFloat
        World units per Angstrom, used to convert Depth Threshold. Molecular Nodes imports structures at 0.1 (1 nm per world unit)

    Inputs
    ------
    i.depth : FloatSocket
        Depth pass from the Render Layers node, in world units
    i.normal : VectorSocket
        Normal pass from the Render Layers node
    i.depth_threshold : FloatSocket
        Depth difference between neighbouring pixels, in Angstrom, above which a line is drawn
    i.use_normals : BooleanSocket
        Also draw lines where the surface normal changes sharply, catching creases between surfaces at the same depth
    i.normal_threshold : FloatSocket
        Change in normal between neighbouring pixels above which a line is drawn, when Use Normals is enabled
    i.size : IntegerSocket
        Line thickness in pixels
    i.world_scale : FloatSocket
        World units per Angstrom, used to convert Depth Threshold. Molecular Nodes imports structures at 0.1 (1 nm per world unit)

    Outputs
    -------
    o.mask : FloatSocket
        Mask
    """

    _name = "Composite Outline Mask"
    _asset_name = "Composite Outline Mask"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "FILTER"
    _tree_properties = {
        "description": "Line mask from the Depth and Normal render passes: a Sobel edge detector on depth, optionally combined with one on normals, anti-aliased and dilated to the requested pixel width"
    }

    class _Inputs(SocketAccessor):
        depth: FloatSocket
        """Depth pass from the Render Layers node, in world units"""
        normal: VectorSocket
        """Normal pass from the Render Layers node"""
        depth_threshold: FloatSocket
        """Depth difference between neighbouring pixels, in Angstrom, above which a line is drawn"""
        use_normals: BooleanSocket
        """Also draw lines where the surface normal changes sharply, catching creases between surfaces at the same depth"""
        normal_threshold: FloatSocket
        """Change in normal between neighbouring pixels above which a line is drawn, when Use Normals is enabled"""
        size: IntegerSocket
        """Line thickness in pixels"""
        world_scale: FloatSocket
        """World units per Angstrom, used to convert Depth Threshold. Molecular Nodes imports structures at 0.1 (1 nm per world unit)"""

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
        normal: InputVector = None,
        depth_threshold: InputFloat = 6.0,
        use_normals: InputBoolean = False,
        normal_threshold: InputFloat = 3.0,
        size: InputInteger = 1,
        world_scale: InputFloat = 0.1,
    ):
        super().__init__(
            **{
                "Depth": depth,
                "Normal": normal,
                "Depth Threshold": depth_threshold,
                "Use Normals": use_normals,
                "Normal Threshold": normal_threshold,
                "Size": size,
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
        normal = tree.inputs.vector(
            "Normal",
            (0.0, 0.0, 0.0),
            description="Normal pass from the Render Layers node",
            hide_value=True,
        )
        depth_threshold = tree.inputs.float(
            "Depth Threshold",
            6.0,
            description="Depth difference between neighbouring pixels, in Angstrom, above which a line is drawn",
            min_value=0.0,
            max_value=10_000.0,
        )
        use_normals = tree.inputs.boolean(
            "Use Normals",
            False,
            description="Also draw lines where the surface normal changes sharply, catching creases between surfaces at the same depth",
        )
        normal_threshold = tree.inputs.float(
            "Normal Threshold",
            3.0,
            description="Change in normal between neighbouring pixels above which a line is drawn, when Use Normals is enabled",
            min_value=0.0,
            max_value=10_000.0,
        )
        size = tree.inputs.integer(
            "Size", 1, description="Line thickness in pixels", min_value=1, max_value=20
        )
        world_scale = tree.inputs.float(
            "World Scale",
            0.1,
            description="World units per Angstrom, used to convert Depth Threshold. Molecular Nodes imports structures at 0.1 (1 nm per world unit)",
            min_value=0.0,
            max_value=10_000.0,
        )
        mask = tree.outputs.float("Mask")

        with c.Frame("Depth edges"):
            math_1 = g.Math.greater_than(
                c.Filter(image=depth.min(10_000.0), type="Sobel"),
                depth_threshold * world_scale,
            )
            _string = g.String(
                string="The depth pass is unbounded where nothing was rendered, so it is clipped before the Sobel filter to keep the gradient finite. Silhouettes against the background still exceed any sensible threshold, so they are outlined."
            )
        with c.Frame("Normal edges"):
            math_2 = g.Math.greater_than(
                g.VectorMath.length(c.Filter(image=normal, type="Sobel")).o.value,
                normal_threshold,
            )
            math_3 = math_1.o.value.max(math_2.o.value * use_normals)
        with c.Frame("Line weight"):
            dilate_erode = c.DilateErode(
                mask=c.AntiAliasing(image=math_3, threshold=0.2),
                size=size - 1.0,
                type="Distance",
            )
            anti_aliasing = c.AntiAliasing(image=dilate_erode, threshold=0.2)

        anti_aliasing >> mask


ASSET = CompositeOutlineMask

ASSET_METADATA = {
    "description": "Line mask from the Depth and Normal render passes",
    "catalog_id": "441e6ca5-e514-4e77-a3cd-25fc1a2e08ae",
}
