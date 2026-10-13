# Node-group asset "Illustrate" (CompositorNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import CompositorNodeTree
from nodebpy import TreeBuilder
from nodebpy import compositor as c
from nodebpy.builder import (
    AssetCompositorGroup,
    BooleanSocket,
    ColorSocket,
    FloatSocket,
    IntegerSocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import (
    InputBoolean,
    InputColor,
    InputFloat,
    InputInteger,
    InputMenu,
    InputVector,
)
from .depth_outline import DepthOutline
from .outline import Outline


class Illustrate(AssetCompositorGroup):
    """
    Illustrative, flat-shaded composite of the render passes: a base colour shaded by the ambient occlusion or shadow pass, outlines from depth and normals, an optional flat background and a depth fade of the alpha

    Parameters
    ----------
    image : InputColor
        Combined pass from the Render Layers node
    alpha : InputFloat
        Alpha pass from the Render Layers node
    depth : InputFloat
        Depth pass from the Render Layers node, in world units
    normal : InputVector
        Normal pass from the Render Layers node
    shading : InputBoolean
        Darken the base colour by the ambient occlusion or shadow pass
    shadow : InputColor
        Shadow pass from the Render Layers node
    shading_intensity : InputFloat
        How dark the fully occluded or shadowed pixels become
    shading_smoothing : InputMenu | Literal["None", "Kuwahara", "Blur"]
        Filter applied to the shading pass before it is used. Kuwahara gives a painterly, flat result; Blur a soft one
    smoothing_size : InputFloat
        Size of the smoothing filter in pixels
    kuwahara_type : InputMenu | Literal["Classic", "Anisotropic"]
        Kuwahara Type
    outline : InputBoolean
        Draw lines where the depth or normal changes sharply
    outline_color : InputColor
        Colour of the lines
    outline_size : InputInteger
        Line thickness in pixels
    outline_method : InputMenu | Literal["Sobel", "Goodsell"]
        Sobel draws hard lines wherever the depth gradient passes Depth Threshold; Goodsell draws soft lines on the nearer side of each depth step, as in Illustrate
    depth_threshold : InputFloat
        Depth gradient, in world units, above which a line is drawn, for the Sobel method
    outline_low : InputFloat
        Weighted sum of the clamped depth steps, in world units, at which a line starts to appear
    outline_high : InputFloat
        Weighted sum of the clamped depth steps, in world units, at which a line is fully opaque
    outline_max_diff : InputFloat
        Depth steps larger than this, in world units, count as this much; a wider range emphasises larger features
    outline_smooth : InputBoolean
        Average the outline over the 3x3 neighbourhood where most pixels carry line signal
    outline_normals : InputBoolean
        Also draw lines where the surface normal changes sharply, catching creases between surfaces at the same depth
    outline_normal_threshold : InputFloat
        Change in normal between neighbouring pixels above which a line is drawn, when Outline Normals is enabled
    background : InputBoolean
        Composite the result over a flat background colour
    background_color : InputColor
        Flat colour behind the structure

    Inputs
    ------
    i.image : ColorSocket
        Combined pass from the Render Layers node
    i.alpha : FloatSocket
        Alpha pass from the Render Layers node
    i.depth : FloatSocket
        Depth pass from the Render Layers node, in world units
    i.normal : VectorSocket
        Normal pass from the Render Layers node
    i.shading : BooleanSocket
        Darken the base colour by the ambient occlusion or shadow pass
    i.shadow : ColorSocket
        Shadow pass from the Render Layers node
    i.shading_intensity : FloatSocket
        How dark the fully occluded or shadowed pixels become
    i.shading_smoothing : MenuSocket
        Filter applied to the shading pass before it is used. Kuwahara gives a painterly, flat result; Blur a soft one
    i.smoothing_size : FloatSocket
        Size of the smoothing filter in pixels
    i.kuwahara_type : MenuSocket
        Kuwahara Type
    i.outline : BooleanSocket
        Draw lines where the depth or normal changes sharply
    i.outline_color : ColorSocket
        Colour of the lines
    i.outline_size : IntegerSocket
        Line thickness in pixels
    i.outline_method : MenuSocket
        Sobel draws hard lines wherever the depth gradient passes Depth Threshold; Goodsell draws soft lines on the nearer side of each depth step, as in Illustrate
    i.depth_threshold : FloatSocket
        Depth gradient, in world units, above which a line is drawn, for the Sobel method
    i.outline_low : FloatSocket
        Weighted sum of the clamped depth steps, in world units, at which a line starts to appear
    i.outline_high : FloatSocket
        Weighted sum of the clamped depth steps, in world units, at which a line is fully opaque
    i.outline_max_diff : FloatSocket
        Depth steps larger than this, in world units, count as this much; a wider range emphasises larger features
    i.outline_smooth : BooleanSocket
        Average the outline over the 3x3 neighbourhood where most pixels carry line signal
    i.outline_normals : BooleanSocket
        Also draw lines where the surface normal changes sharply, catching creases between surfaces at the same depth
    i.outline_normal_threshold : FloatSocket
        Change in normal between neighbouring pixels above which a line is drawn, when Outline Normals is enabled
    i.background : BooleanSocket
        Composite the result over a flat background colour
    i.background_color : ColorSocket
        Flat colour behind the structure

    Outputs
    -------
    o.image : ColorSocket
        Image
    """

    _name = "Illustrate"
    _asset_name = "Illustrate"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "FILTER"
    _tree_properties = {
        "description": "Illustrative, flat-shaded composite of the render passes: a base colour shaded by the ambient occlusion or shadow pass, outlines from depth and normals, an optional flat background and a depth fade of the alpha"
    }

    class _Inputs(SocketAccessor):
        image: ColorSocket
        """Combined pass from the Render Layers node"""
        alpha: FloatSocket
        """Alpha pass from the Render Layers node"""
        depth: FloatSocket
        """Depth pass from the Render Layers node, in world units"""
        normal: VectorSocket
        """Normal pass from the Render Layers node"""
        shading: BooleanSocket
        """Darken the base colour by the ambient occlusion or shadow pass"""
        shadow: ColorSocket
        """Shadow pass from the Render Layers node"""
        shading_intensity: FloatSocket
        """How dark the fully occluded or shadowed pixels become"""
        shading_smoothing: MenuSocket
        """Filter applied to the shading pass before it is used. Kuwahara gives a painterly, flat result; Blur a soft one"""
        smoothing_size: FloatSocket
        """Size of the smoothing filter in pixels"""
        kuwahara_type: MenuSocket
        """Kuwahara Type"""
        outline: BooleanSocket
        """Draw lines where the depth or normal changes sharply"""
        outline_color: ColorSocket
        """Colour of the lines"""
        outline_size: IntegerSocket
        """Line thickness in pixels"""
        outline_method: MenuSocket
        """Sobel draws hard lines wherever the depth gradient passes Depth Threshold; Goodsell draws soft lines on the nearer side of each depth step, as in Illustrate"""
        depth_threshold: FloatSocket
        """Depth gradient, in world units, above which a line is drawn, for the Sobel method"""
        outline_low: FloatSocket
        """Weighted sum of the clamped depth steps, in world units, at which a line starts to appear"""
        outline_high: FloatSocket
        """Weighted sum of the clamped depth steps, in world units, at which a line is fully opaque"""
        outline_max_diff: FloatSocket
        """Depth steps larger than this, in world units, count as this much; a wider range emphasises larger features"""
        outline_smooth: BooleanSocket
        """Average the outline over the 3x3 neighbourhood where most pixels carry line signal"""
        outline_normals: BooleanSocket
        """Also draw lines where the surface normal changes sharply, catching creases between surfaces at the same depth"""
        outline_normal_threshold: FloatSocket
        """Change in normal between neighbouring pixels above which a line is drawn, when Outline Normals is enabled"""
        background: BooleanSocket
        """Composite the result over a flat background colour"""
        background_color: ColorSocket
        """Flat colour behind the structure"""

    class _Outputs(SocketAccessor):
        image: ColorSocket
        """Image"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        image: InputColor = None,
        alpha: InputFloat = 1.0,
        depth: InputFloat = 0.0,
        normal: InputVector = None,
        shading: InputBoolean = True,
        shadow: InputColor = None,
        shading_intensity: InputFloat = 1.0,
        shading_smoothing: InputMenu | Literal["None", "Kuwahara", "Blur"] = "Kuwahara",
        smoothing_size: InputFloat = 10.0,
        kuwahara_type: InputMenu | Literal["Classic", "Anisotropic"] = "Classic",
        outline: InputBoolean = True,
        outline_color: InputColor = None,
        outline_size: InputInteger = 1,
        outline_method: InputMenu | Literal["Sobel", "Goodsell"] = "Sobel",
        depth_threshold: InputFloat = 1.0,
        outline_low: InputFloat = 0.3,
        outline_high: InputFloat = 1.0,
        outline_max_diff: InputFloat = 0.5,
        outline_smooth: InputBoolean = True,
        outline_normals: InputBoolean = False,
        outline_normal_threshold: InputFloat = 3.0,
        background: InputBoolean = False,
        background_color: InputColor = None,
    ):
        super().__init__(
            **{
                "Image": image,
                "Alpha": alpha,
                "Depth": depth,
                "Normal": normal,
                "Shading": shading,
                "Shadow": shadow,
                "Shading Intensity": shading_intensity,
                "Shading Smoothing": shading_smoothing,
                "Smoothing Size": smoothing_size,
                "Kuwahara Type": kuwahara_type,
                "Outline": outline,
                "Outline Color": outline_color,
                "Outline Size": outline_size,
                "Outline Method": outline_method,
                "Depth Threshold": depth_threshold,
                "Outline Low": outline_low,
                "Outline High": outline_high,
                "Outline Max Diff": outline_max_diff,
                "Outline Smooth": outline_smooth,
                "Outline Normals": outline_normals,
                "Outline Normal Threshold": outline_normal_threshold,
                "Background": background,
                "Background Color": background_color,
            }
        )

    def _build_group(self, tree: TreeBuilder[CompositorNodeTree]) -> None:
        image = tree.inputs.color(
            "Image",
            (0.0, 0.0, 0.0, 1.0),
            description="Combined pass from the Render Layers node",
            hide_value=True,
        )
        alpha = tree.inputs.float(
            "Alpha",
            1.0,
            description="Alpha pass from the Render Layers node",
            min_value=0.0,
            max_value=1.0,
            hide_value=True,
            subtype="FACTOR",
        )
        depth = tree.inputs.float(
            "Depth",
            0.0,
            description="Depth pass from the Render Layers node, in world units",
            hide_value=True,
            subtype="DISTANCE",
        )
        normal = tree.inputs.vector(
            "Normal",
            (0.0, 0.0, 0.0),
            description="Normal pass from the Render Layers node",
            hide_value=True,
        )
        with tree.inputs.panel("Shading"):
            shading = tree.inputs.boolean(
                "Shading",
                True,
                description="Darken the base colour by the ambient occlusion or shadow pass",
                is_panel_toggle=True,
            )
            shadow = tree.inputs.color(
                "Shadow",
                (1.0, 1.0, 1.0, 1.0),
                description="Shadow pass from the Render Layers node",
                hide_value=True,
            )
            shading_intensity = tree.inputs.float(
                "Shading Intensity",
                1.0,
                description="How dark the fully occluded or shadowed pixels become",
                min_value=0.0,
                max_value=1.0,
                subtype="FACTOR",
            )
            shading_smoothing = tree.inputs.menu(
                "Shading Smoothing",
                "Kuwahara",
                description="Filter applied to the shading pass before it is used. Kuwahara gives a painterly, flat result; Blur a soft one",
                optional_label=True,
            )
            smoothing_size = tree.inputs.float(
                "Smoothing Size",
                10.0,
                description="Size of the smoothing filter in pixels",
                min_value=0.0,
                max_value=100.0,
            )
            kuwahara_type = tree.inputs.menu(
                "Kuwahara Type", "Classic", optional_label=True
            )
        with tree.inputs.panel("Outline"):
            outline = tree.inputs.boolean(
                "Outline",
                True,
                description="Draw lines where the depth or normal changes sharply",
                is_panel_toggle=True,
            )
            outline_color = tree.inputs.color(
                "Outline Color", (0.0, 0.0, 0.0, 1.0), description="Colour of the lines"
            )
            outline_size = tree.inputs.integer(
                "Outline Size",
                1,
                description="Line thickness in pixels",
                min_value=1,
                max_value=20,
            )
            outline_method = tree.inputs.menu(
                "Outline Method",
                "Sobel",
                description="Sobel draws hard lines wherever the depth gradient passes Depth Threshold; Goodsell draws soft lines on the nearer side of each depth step, as in Illustrate",
            )
            depth_threshold = tree.inputs.float(
                "Depth Threshold",
                1.0,
                description="Depth gradient, in world units, above which a line is drawn, for the Sobel method",
                min_value=0.0,
                max_value=10_000.0,
            )
            outline_low = tree.inputs.float(
                "Outline Low",
                0.3,
                description="Weighted sum of the clamped depth steps, in world units, at which a line starts to appear",
                min_value=0.0,
                max_value=10_000.0,
            )
            outline_high = tree.inputs.float(
                "Outline High",
                1.0,
                description="Weighted sum of the clamped depth steps, in world units, at which a line is fully opaque",
                min_value=0.0,
                max_value=10_000.0,
            )
            outline_max_diff = tree.inputs.float(
                "Outline Max Diff",
                0.5,
                description="Depth steps larger than this, in world units, count as this much; a wider range emphasises larger features",
                min_value=0.0,
                max_value=10_000.0,
            )
            outline_smooth = tree.inputs.boolean(
                "Outline Smooth",
                True,
                description="Average the outline over the 3x3 neighbourhood where most pixels carry line signal",
            )
            outline_normals = tree.inputs.boolean(
                "Outline Normals",
                False,
                description="Also draw lines where the surface normal changes sharply, catching creases between surfaces at the same depth",
            )
            outline_normal_threshold = tree.inputs.float(
                "Outline Normal Threshold",
                3.0,
                description="Change in normal between neighbouring pixels above which a line is drawn, when Outline Normals is enabled",
                min_value=0.0,
                max_value=10_000.0,
            )
        with tree.inputs.panel("Background"):
            background = tree.inputs.boolean(
                "Background",
                False,
                description="Composite the result over a flat background colour",
                is_panel_toggle=True,
            )
            background_color = tree.inputs.color(
                "Background Color",
                (1.0, 1.0, 1.0, 1.0),
                description="Flat colour behind the structure",
            )
        image_1 = tree.outputs.color("Image", (1.0, 1.0, 1.0, 1.0))

        with c.Frame("Shading"):
            kuwahara = c.Kuwahara(
                image=shadow,
                size=smoothing_size,
                type=kuwahara_type,
                uniformity=5,
                sharpness=0.5,
            )
            menu_switch = c.MenuSwitch.color(
                shading_smoothing,
                {
                    "None": shadow,
                    "Kuwahara": kuwahara,
                    "Blur": c.Blur(image=shadow, size=smoothing_size),
                },
            )
            mix = (shading_intensity * c.InvertColor(menu_switch.o.output)).mix.color(
                image, (0.0, 0.0, 0.0, 1.0)
            )
            _string = c.String(
                "The chosen pass is smoothed, inverted to an occlusion amount, scaled by the intensity and used to mix the base colour towards black."
            )
            set_alpha = c.SetAlpha(image=shading.switch.color(image, mix), alpha=alpha)
        with c.Frame("Outline"):
            depth_outline = DepthOutline(
                depth=depth,
                method=outline_method,
                size=outline_size,
                threshold=depth_threshold,
                low=outline_low,
                high=outline_high,
                max_diff=outline_max_diff,
                smooth=outline_smooth,
                normals=outline_normals,
                normals_pass=normal,
                normals_threshold=outline_normal_threshold,
            )
            _string_1 = c.String(
                "The alpha is replaced once here: the Diffuse Color pass has no alpha and mixing keeps only the first input's, so the render alpha (or its depth fade) is restored before the lines are drawn on top."
            )
            switch = outline.switch.color(
                set_alpha,
                Outline(image=set_alpha, mask=depth_outline, color=outline_color),
            )
        with c.Frame("Background"):
            (
                background.switch.color(
                    switch,
                    c.AlphaOver(
                        background=background_color,
                        foreground=switch,
                        straight_alpha=True,
                    ),
                )
                >> image_1
            )


ASSET = Illustrate

ASSET_METADATA = {
    "description": "Illustrative, flat-shaded composite of the render passes with outlines",
    "catalog_id": "874a72a6-7cf2-456a-ab07-e79f5a1e4294",
}
