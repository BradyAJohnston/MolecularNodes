# Node-group asset "Composite Illustrative" (CompositorNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import CompositorNodeTree
from nodebpy import TreeBuilder
from nodebpy import compositor as c
from nodebpy import geometry as g
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
from .composite_outline import CompositeOutline
from .composite_outline_mask import CompositeOutlineMask


class CompositeIllustrative(AssetCompositorGroup):
    """
    Illustrative, flat-shaded composite of the render passes: a base colour shaded by the ambient occlusion or shadow pass, outlines from depth and normals, an optional flat background and a depth fade of the alpha

    Parameters
    ----------
    base_color : InputMenu | Literal["Image", "Diffuse Color"]
        Whether to shade the rendered image or the flat Diffuse Color pass, which ignores lighting and materials entirely
    image : InputColor
        Combined pass from the Render Layers node
    alpha : InputFloat
        Alpha pass from the Render Layers node
    depth : InputFloat
        Depth pass from the Render Layers node, in world units
    normal : InputVector
        Normal pass from the Render Layers node
    diffuse_color : InputColor
        Diffuse Color pass from the Render Layers node
    ambient_occlusion : InputColor
        Ambient Occlusion pass from the Render Layers node
    shadow : InputColor
        Shadow pass from the Render Layers node
    world_scale : InputFloat
        World units per Angstrom, used to convert the Angstrom inputs. Molecular Nodes imports structures at 0.1 (1 nm per world unit)
    shading : InputBoolean
        Darken the base colour by the ambient occlusion or shadow pass
    shading_source : InputMenu | Literal["Ambient Occlusion", "Shadow"]
        Which pass darkens the base colour
    shading_intensity : InputFloat
        How dark the fully occluded or shadowed pixels become
    shading_smoothing : InputMenu | Literal["None", "Kuwahara", "Blur"]
        Filter applied to the shading pass before it is used. Kuwahara gives a painterly, flat result; Blur a soft one
    smoothing_size : InputFloat
        Size of the smoothing filter in pixels
    outline : InputBoolean
        Draw lines where the depth or normal changes sharply
    outline_color : InputColor
        Colour of the lines
    outline_size : InputInteger
        Line thickness in pixels
    depth_threshold : InputFloat
        Depth difference between neighbouring pixels, in Angstrom, above which a line is drawn
    use_normals : InputBoolean
        Also draw lines where the surface normal changes sharply, catching creases between surfaces at the same depth
    normal_threshold : InputFloat
        Change in normal between neighbouring pixels above which a line is drawn, when Use Normals is enabled
    background : InputBoolean
        Composite the result over a flat background colour
    background_color : InputColor
        Flat colour behind the structure
    depth_fade : InputBoolean
        Fade the alpha to zero with distance from the camera
    alpha_distance : InputFloat
        Distance from the camera, in Angstrom, at the middle of the fade
    alpha_falloff : InputFloat
        Half-width of the fade in Angstrom: fully opaque this far in front of Alpha Distance, fully transparent this far behind it
    alpha_power : InputFloat
        Exponent applied to the fade; lower values keep more of the structure visible before it drops away

    Inputs
    ------
    i.base_color : MenuSocket
        Whether to shade the rendered image or the flat Diffuse Color pass, which ignores lighting and materials entirely
    i.image : ColorSocket
        Combined pass from the Render Layers node
    i.alpha : FloatSocket
        Alpha pass from the Render Layers node
    i.depth : FloatSocket
        Depth pass from the Render Layers node, in world units
    i.normal : VectorSocket
        Normal pass from the Render Layers node
    i.diffuse_color : ColorSocket
        Diffuse Color pass from the Render Layers node
    i.ambient_occlusion : ColorSocket
        Ambient Occlusion pass from the Render Layers node
    i.shadow : ColorSocket
        Shadow pass from the Render Layers node
    i.world_scale : FloatSocket
        World units per Angstrom, used to convert the Angstrom inputs. Molecular Nodes imports structures at 0.1 (1 nm per world unit)
    i.shading : BooleanSocket
        Darken the base colour by the ambient occlusion or shadow pass
    i.shading_source : MenuSocket
        Which pass darkens the base colour
    i.shading_intensity : FloatSocket
        How dark the fully occluded or shadowed pixels become
    i.shading_smoothing : MenuSocket
        Filter applied to the shading pass before it is used. Kuwahara gives a painterly, flat result; Blur a soft one
    i.smoothing_size : FloatSocket
        Size of the smoothing filter in pixels
    i.outline : BooleanSocket
        Draw lines where the depth or normal changes sharply
    i.outline_color : ColorSocket
        Colour of the lines
    i.outline_size : IntegerSocket
        Line thickness in pixels
    i.depth_threshold : FloatSocket
        Depth difference between neighbouring pixels, in Angstrom, above which a line is drawn
    i.use_normals : BooleanSocket
        Also draw lines where the surface normal changes sharply, catching creases between surfaces at the same depth
    i.normal_threshold : FloatSocket
        Change in normal between neighbouring pixels above which a line is drawn, when Use Normals is enabled
    i.background : BooleanSocket
        Composite the result over a flat background colour
    i.background_color : ColorSocket
        Flat colour behind the structure
    i.depth_fade : BooleanSocket
        Fade the alpha to zero with distance from the camera
    i.alpha_distance : FloatSocket
        Distance from the camera, in Angstrom, at the middle of the fade
    i.alpha_falloff : FloatSocket
        Half-width of the fade in Angstrom: fully opaque this far in front of Alpha Distance, fully transparent this far behind it
    i.alpha_power : FloatSocket
        Exponent applied to the fade; lower values keep more of the structure visible before it drops away

    Outputs
    -------
    o.image : ColorSocket
        Image
    """

    _name = "Composite Illustrative"
    _asset_name = "Composite Illustrative"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "FILTER"
    _tree_properties = {
        "description": "Illustrative, flat-shaded composite of the render passes: a base colour shaded by the ambient occlusion or shadow pass, outlines from depth and normals, an optional flat background and a depth fade of the alpha"
    }

    class _Inputs(SocketAccessor):
        base_color: MenuSocket
        """Whether to shade the rendered image or the flat Diffuse Color pass, which ignores lighting and materials entirely"""
        image: ColorSocket
        """Combined pass from the Render Layers node"""
        alpha: FloatSocket
        """Alpha pass from the Render Layers node"""
        depth: FloatSocket
        """Depth pass from the Render Layers node, in world units"""
        normal: VectorSocket
        """Normal pass from the Render Layers node"""
        diffuse_color: ColorSocket
        """Diffuse Color pass from the Render Layers node"""
        ambient_occlusion: ColorSocket
        """Ambient Occlusion pass from the Render Layers node"""
        shadow: ColorSocket
        """Shadow pass from the Render Layers node"""
        world_scale: FloatSocket
        """World units per Angstrom, used to convert the Angstrom inputs. Molecular Nodes imports structures at 0.1 (1 nm per world unit)"""
        shading: BooleanSocket
        """Darken the base colour by the ambient occlusion or shadow pass"""
        shading_source: MenuSocket
        """Which pass darkens the base colour"""
        shading_intensity: FloatSocket
        """How dark the fully occluded or shadowed pixels become"""
        shading_smoothing: MenuSocket
        """Filter applied to the shading pass before it is used. Kuwahara gives a painterly, flat result; Blur a soft one"""
        smoothing_size: FloatSocket
        """Size of the smoothing filter in pixels"""
        outline: BooleanSocket
        """Draw lines where the depth or normal changes sharply"""
        outline_color: ColorSocket
        """Colour of the lines"""
        outline_size: IntegerSocket
        """Line thickness in pixels"""
        depth_threshold: FloatSocket
        """Depth difference between neighbouring pixels, in Angstrom, above which a line is drawn"""
        use_normals: BooleanSocket
        """Also draw lines where the surface normal changes sharply, catching creases between surfaces at the same depth"""
        normal_threshold: FloatSocket
        """Change in normal between neighbouring pixels above which a line is drawn, when Use Normals is enabled"""
        background: BooleanSocket
        """Composite the result over a flat background colour"""
        background_color: ColorSocket
        """Flat colour behind the structure"""
        depth_fade: BooleanSocket
        """Fade the alpha to zero with distance from the camera"""
        alpha_distance: FloatSocket
        """Distance from the camera, in Angstrom, at the middle of the fade"""
        alpha_falloff: FloatSocket
        """Half-width of the fade in Angstrom: fully opaque this far in front of Alpha Distance, fully transparent this far behind it"""
        alpha_power: FloatSocket
        """Exponent applied to the fade; lower values keep more of the structure visible before it drops away"""

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
        base_color: InputMenu | Literal["Image", "Diffuse Color"] = "Image",
        image: InputColor = None,
        alpha: InputFloat = 1.0,
        depth: InputFloat = 0.0,
        normal: InputVector = None,
        diffuse_color: InputColor = None,
        ambient_occlusion: InputColor = None,
        shadow: InputColor = None,
        world_scale: InputFloat = 0.1,
        shading: InputBoolean = True,
        shading_source: InputMenu
        | Literal["Ambient Occlusion", "Shadow"] = "Ambient Occlusion",
        shading_intensity: InputFloat = 1.0,
        shading_smoothing: InputMenu | Literal["None", "Kuwahara", "Blur"] = "Kuwahara",
        smoothing_size: InputFloat = 5.0,
        outline: InputBoolean = True,
        outline_color: InputColor = None,
        outline_size: InputInteger = 2,
        depth_threshold: InputFloat = 6.0,
        use_normals: InputBoolean = False,
        normal_threshold: InputFloat = 3.0,
        background: InputBoolean = False,
        background_color: InputColor = None,
        depth_fade: InputBoolean = False,
        alpha_distance: InputFloat = 100.0,
        alpha_falloff: InputFloat = 20.0,
        alpha_power: InputFloat = 0.25,
    ):
        super().__init__(
            **{
                "Base Color": base_color,
                "Image": image,
                "Alpha": alpha,
                "Depth": depth,
                "Normal": normal,
                "Diffuse Color": diffuse_color,
                "Ambient Occlusion": ambient_occlusion,
                "Shadow": shadow,
                "World Scale": world_scale,
                "Shading": shading,
                "Shading Source": shading_source,
                "Shading Intensity": shading_intensity,
                "Shading Smoothing": shading_smoothing,
                "Smoothing Size": smoothing_size,
                "Outline": outline,
                "Outline Color": outline_color,
                "Outline Size": outline_size,
                "Depth Threshold": depth_threshold,
                "Use Normals": use_normals,
                "Normal Threshold": normal_threshold,
                "Background": background,
                "Background Color": background_color,
                "Depth Fade": depth_fade,
                "Alpha Distance": alpha_distance,
                "Alpha Falloff": alpha_falloff,
                "Alpha Power": alpha_power,
            }
        )

    def _build_group(self, tree: TreeBuilder[CompositorNodeTree]) -> None:
        base_color = tree.inputs.menu(
            "Base Color",
            description="Whether to shade the rendered image or the flat Diffuse Color pass, which ignores lighting and materials entirely",
            expanded=True,
        )
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
        )
        normal = tree.inputs.vector(
            "Normal",
            (0.0, 0.0, 0.0),
            description="Normal pass from the Render Layers node",
            hide_value=True,
        )
        diffuse_color = tree.inputs.color(
            "Diffuse Color",
            (0.8, 0.8, 0.8, 1.0),
            description="Diffuse Color pass from the Render Layers node",
            hide_value=True,
        )
        ambient_occlusion = tree.inputs.color(
            "Ambient Occlusion",
            (1.0, 1.0, 1.0, 1.0),
            description="Ambient Occlusion pass from the Render Layers node",
            hide_value=True,
        )
        shadow = tree.inputs.color(
            "Shadow",
            (1.0, 1.0, 1.0, 1.0),
            description="Shadow pass from the Render Layers node",
            hide_value=True,
        )
        world_scale = tree.inputs.float(
            "World Scale",
            0.1,
            description="World units per Angstrom, used to convert the Angstrom inputs. Molecular Nodes imports structures at 0.1 (1 nm per world unit)",
            min_value=0.0,
            max_value=10_000.0,
        )
        with tree.inputs.panel("Shading"):
            shading = tree.inputs.boolean(
                "Shading",
                True,
                description="Darken the base colour by the ambient occlusion or shadow pass",
                is_panel_toggle=True,
            )
            shading_source = tree.inputs.menu(
                "Shading Source",
                description="Which pass darkens the base colour",
                expanded=True,
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
                description="Filter applied to the shading pass before it is used. Kuwahara gives a painterly, flat result; Blur a soft one",
                expanded=True,
            )
            smoothing_size = tree.inputs.float(
                "Smoothing Size",
                5.0,
                description="Size of the smoothing filter in pixels",
                min_value=0.0,
                max_value=100.0,
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
                2,
                description="Line thickness in pixels",
                min_value=1,
                max_value=20,
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
        with tree.inputs.panel("Depth Fade", default_closed=True):
            depth_fade = tree.inputs.boolean(
                "Depth Fade",
                False,
                description="Fade the alpha to zero with distance from the camera",
                is_panel_toggle=True,
            )
            alpha_distance = tree.inputs.float(
                "Alpha Distance",
                100.0,
                description="Distance from the camera, in Angstrom, at the middle of the fade",
                min_value=0.0,
                max_value=100_000.0,
            )
            alpha_falloff = tree.inputs.float(
                "Alpha Falloff",
                20.0,
                description="Half-width of the fade in Angstrom: fully opaque this far in front of Alpha Distance, fully transparent this far behind it",
                min_value=0.0,
                max_value=100_000.0,
            )
            alpha_power = tree.inputs.float(
                "Alpha Power",
                0.25,
                description="Exponent applied to the fade; lower values keep more of the structure visible before it drops away",
                min_value=0.0,
                max_value=1.0,
                subtype="FACTOR",
            )
        image_1 = tree.outputs.color("Image", (1.0, 1.0, 1.0, 1.0))

        with c.Frame("Depth fade"):
            math_1 = alpha_distance * world_scale
            math_2 = alpha_falloff.max(0.001) * world_scale
            math_3 = depth.map_range(
                math_1 + math_2, math_1 - math_2, interpolation_type="SMOOTHERSTEP"
            ) ** alpha_power.max(0.001)
            _string = g.String(
                string="Alpha is kept in front of Alpha Distance minus Alpha Falloff and fades to zero by Alpha Distance plus Alpha Falloff with a smootherstep ramp raised to Alpha Power. It never exceeds the render's own alpha."
            )
            switch = c.Switch(switch=depth_fade, off=alpha, on=math_3.min(alpha))
        with c.Frame("Base colour"):
            menu_switch = c.MenuSwitch.color(
                base_color, {"Image": image, "Diffuse Color": diffuse_color}
            )
        with c.Frame("Shading"):
            menu_switch_1 = c.MenuSwitch.color(
                shading_source,
                {"Ambient Occlusion": ambient_occlusion, "Shadow": shadow},
            )
            menu_switch_2 = c.MenuSwitch.color(
                shading_smoothing,
                {
                    "None": menu_switch_1.o.output,
                    "Kuwahara": c.Kuwahara(
                        image=menu_switch_1.o.output, size=smoothing_size, sharpness=0.5
                    ),
                    "Blur": c.Blur(image=menu_switch_1.o.output, size=smoothing_size),
                },
            )
            mix = (
                shading_intensity * c.InvertColor(color=menu_switch_2.o.output)
            ).mix.color(menu_switch.o.output, (0.0, 0.0, 0.0, 1.0))
            _string_1 = g.String(
                string="The chosen pass is smoothed, inverted to an occlusion amount, scaled by the intensity and used to mix the base colour towards black."
            )
            switch_1 = c.Switch(switch=shading, off=menu_switch.o.output, on=mix)
        with c.Frame("Outline"):
            composite_outline_mask = CompositeOutlineMask(
                depth=depth,
                normal=normal,
                depth_threshold=depth_threshold,
                use_normals=use_normals,
                normal_threshold=normal_threshold,
                size=outline_size,
                world_scale=world_scale,
            )
            _string_2 = g.String(
                string="The alpha is replaced once here: the Diffuse Color pass has no alpha and mixing keeps only the first input's, so the render alpha (or its depth fade) is restored before the lines are drawn on top."
            )
            set_alpha = c.SetAlpha(image=switch_1, alpha=switch, type="Replace Alpha")
            switch_2 = c.Switch(
                switch=outline,
                off=set_alpha,
                on=CompositeOutline(
                    image=set_alpha,
                    mask=composite_outline_mask,
                    line_color=outline_color,
                ),
            )
        with c.Frame("Background"):
            switch_3 = c.Switch(
                switch=background,
                off=switch_2,
                on=c.AlphaOver(
                    background=background_color,
                    foreground=switch_2,
                    straight_alpha=True,
                ),
            )

        switch_3 >> image_1

        base_color.default_value = "Image"
        shading_source.default_value = "Ambient Occlusion"
        shading_smoothing.default_value = "Kuwahara"


ASSET = CompositeIllustrative

ASSET_METADATA = {
    "description": "Illustrative, flat-shaded composite of the render passes with outlines",
    "catalog_id": "441e6ca5-e514-4e77-a3cd-25fc1a2e08ae",
}
