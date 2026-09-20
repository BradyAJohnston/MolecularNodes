# Node-group asset "Composite Illustrate" (CompositorNodeTree), dumped by nodebpy.assets.dump_library.
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
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputColor, InputFloat, InputMenu
from .composite_cone_shadow import CompositeConeShadow
from .composite_contour_outline import CompositeContourOutline
from .composite_depth_fog import CompositeDepthFog
from .composite_id_outline import CompositeIDOutline


class CompositeIllustrate(AssetCompositorGroup):
    """
    David Goodsell's Illustrate composition from the render passes: a flat base colour times the conical soft shadow, blended with depth fog, darkened by the contour, residue and chain outlines, with the outlines extending the alpha. Meant for the Flat material, which has no lighting of its own

    Parameters
    ----------
    base_color : InputMenu | Literal["Diffuse Color", "Image"]
        The flat Diffuse Color pass, which ignores lighting and materials, or the rendered image
    image : InputColor
        Combined pass from the Render Layers node
    alpha : InputFloat
        Alpha pass from the Render Layers node
    depth : InputFloat
        Depth pass from the Render Layers node, in world units
    diffuse_color : InputColor
        Diffuse Color pass from the Render Layers node
    world_scale : InputFloat
        World units per Angstrom, used to convert the Angstrom inputs. Molecular Nodes imports structures at 0.1 (1 nm per world unit)
    shadow : InputBoolean
        Darken each pixel by the surrounding pixels that are closer to the camera (Composite Cone Shadow)
    radius : InputFloat
        Radius of the shadow neighbourhood in pixels
    cone_angle : InputFloat
        Tightness of the shadow cone; higher values keep shadows local to crevices
    min_gap : InputFloat
        Smallest depth gap, in Angstrom, that casts a shadow
    strength : InputFloat
        Darkening when every sample occludes the pixel
    max_darkening : InputFloat
        The shadow factor never drops below this value
    pixel_size : InputFloat
        Width of one pixel in Angstrom at the structure. canvas.compositor.illustrate() sets it from the camera
    fog : InputBoolean
        Blend towards Fog Color with depth (Composite Depth Fog)
    near : InputFloat
        Depth, in world units, of the front of the structure. canvas.compositor.illustrate() sets it from the scene bounds
    far : InputFloat
        Depth, in world units, of the back of the structure
    front_fog : InputFloat
        Fog factor at Near: 1 keeps the colour, 0 is fully fogged
    back_fog : InputFloat
        Fog factor at Far: 1 keeps the colour, 0 is fully fogged
    fog_color : InputColor
        Colour the structure fades towards
    contour_outline : InputBoolean
        Lines on depth steps from a Laplacian of the depth pass (Composite Contour Outline)
    contour_high : InputFloat
        Laplacian value, in Angstrom, at which the line is fully opaque
    smooth : InputBoolean
        Average the line opacity where most of the 3x3 neighbourhood carries signal
    chain_outline : InputBoolean
        Lines between chains from a chain_id AOV pass (Composite ID Outline)
    chain_id : InputFloat
        chain_id AOV pass from the Render Layers node
    chain_low : InputFloat
        Number of neighbours in another chain, of 24, at which the line starts to appear
    chain_high : InputFloat
        Number of neighbours in another chain, of 24, at which the line is fully opaque
    residue_outline : InputBoolean
        Lines between residues from a res_id AOV pass (Composite ID Outline)
    residue_id : InputFloat
        res_id AOV pass from the Render Layers node
    residue_low : InputFloat
        Number of neighbours in another residue, of 24, at which the line starts to appear
    residue_high : InputFloat
        Number of neighbours in another residue, of 24, at which the line is fully opaque
    residue_difference : InputFloat
        Difference in residue number above which two pixels are in different residues; raise it to outline groups of residues

    Inputs
    ------
    i.base_color : MenuSocket
        The flat Diffuse Color pass, which ignores lighting and materials, or the rendered image
    i.image : ColorSocket
        Combined pass from the Render Layers node
    i.alpha : FloatSocket
        Alpha pass from the Render Layers node
    i.depth : FloatSocket
        Depth pass from the Render Layers node, in world units
    i.diffuse_color : ColorSocket
        Diffuse Color pass from the Render Layers node
    i.world_scale : FloatSocket
        World units per Angstrom, used to convert the Angstrom inputs. Molecular Nodes imports structures at 0.1 (1 nm per world unit)
    i.shadow : BooleanSocket
        Darken each pixel by the surrounding pixels that are closer to the camera (Composite Cone Shadow)
    i.radius : FloatSocket
        Radius of the shadow neighbourhood in pixels
    i.cone_angle : FloatSocket
        Tightness of the shadow cone; higher values keep shadows local to crevices
    i.min_gap : FloatSocket
        Smallest depth gap, in Angstrom, that casts a shadow
    i.strength : FloatSocket
        Darkening when every sample occludes the pixel
    i.max_darkening : FloatSocket
        The shadow factor never drops below this value
    i.pixel_size : FloatSocket
        Width of one pixel in Angstrom at the structure. canvas.compositor.illustrate() sets it from the camera
    i.fog : BooleanSocket
        Blend towards Fog Color with depth (Composite Depth Fog)
    i.near : FloatSocket
        Depth, in world units, of the front of the structure. canvas.compositor.illustrate() sets it from the scene bounds
    i.far : FloatSocket
        Depth, in world units, of the back of the structure
    i.front_fog : FloatSocket
        Fog factor at Near: 1 keeps the colour, 0 is fully fogged
    i.back_fog : FloatSocket
        Fog factor at Far: 1 keeps the colour, 0 is fully fogged
    i.fog_color : ColorSocket
        Colour the structure fades towards
    i.contour_outline : BooleanSocket
        Lines on depth steps from a Laplacian of the depth pass (Composite Contour Outline)
    i.contour_high : FloatSocket
        Laplacian value, in Angstrom, at which the line is fully opaque
    i.smooth : BooleanSocket
        Average the line opacity where most of the 3x3 neighbourhood carries signal
    i.chain_outline : BooleanSocket
        Lines between chains from a chain_id AOV pass (Composite ID Outline)
    i.chain_id : FloatSocket
        chain_id AOV pass from the Render Layers node
    i.chain_low : FloatSocket
        Number of neighbours in another chain, of 24, at which the line starts to appear
    i.chain_high : FloatSocket
        Number of neighbours in another chain, of 24, at which the line is fully opaque
    i.residue_outline : BooleanSocket
        Lines between residues from a res_id AOV pass (Composite ID Outline)
    i.residue_id : FloatSocket
        res_id AOV pass from the Render Layers node
    i.residue_low : FloatSocket
        Number of neighbours in another residue, of 24, at which the line starts to appear
    i.residue_high : FloatSocket
        Number of neighbours in another residue, of 24, at which the line is fully opaque
    i.residue_difference : FloatSocket
        Difference in residue number above which two pixels are in different residues; raise it to outline groups of residues

    Outputs
    -------
    o.image : ColorSocket
        Image
    """

    _name = "Composite Illustrate"
    _asset_name = "Composite Illustrate"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "FILTER"
    _tree_properties = {
        "description": "David Goodsell's Illustrate composition from the render passes: a flat base colour times the conical soft shadow, blended with depth fog, darkened by the contour, residue and chain outlines, with the outlines extending the alpha. Meant for the Flat material, which has no lighting of its own"
    }

    class _Inputs(SocketAccessor):
        base_color: MenuSocket
        """The flat Diffuse Color pass, which ignores lighting and materials, or the rendered image"""
        image: ColorSocket
        """Combined pass from the Render Layers node"""
        alpha: FloatSocket
        """Alpha pass from the Render Layers node"""
        depth: FloatSocket
        """Depth pass from the Render Layers node, in world units"""
        diffuse_color: ColorSocket
        """Diffuse Color pass from the Render Layers node"""
        world_scale: FloatSocket
        """World units per Angstrom, used to convert the Angstrom inputs. Molecular Nodes imports structures at 0.1 (1 nm per world unit)"""
        shadow: BooleanSocket
        """Darken each pixel by the surrounding pixels that are closer to the camera (Composite Cone Shadow)"""
        radius: FloatSocket
        """Radius of the shadow neighbourhood in pixels"""
        cone_angle: FloatSocket
        """Tightness of the shadow cone; higher values keep shadows local to crevices"""
        min_gap: FloatSocket
        """Smallest depth gap, in Angstrom, that casts a shadow"""
        strength: FloatSocket
        """Darkening when every sample occludes the pixel"""
        max_darkening: FloatSocket
        """The shadow factor never drops below this value"""
        pixel_size: FloatSocket
        """Width of one pixel in Angstrom at the structure. canvas.compositor.illustrate() sets it from the camera"""
        fog: BooleanSocket
        """Blend towards Fog Color with depth (Composite Depth Fog)"""
        near: FloatSocket
        """Depth, in world units, of the front of the structure. canvas.compositor.illustrate() sets it from the scene bounds"""
        far: FloatSocket
        """Depth, in world units, of the back of the structure"""
        front_fog: FloatSocket
        """Fog factor at Near: 1 keeps the colour, 0 is fully fogged"""
        back_fog: FloatSocket
        """Fog factor at Far: 1 keeps the colour, 0 is fully fogged"""
        fog_color: ColorSocket
        """Colour the structure fades towards"""
        contour_outline: BooleanSocket
        """Lines on depth steps from a Laplacian of the depth pass (Composite Contour Outline)"""
        contour_high: FloatSocket
        """Laplacian value, in Angstrom, at which the line is fully opaque"""
        smooth: BooleanSocket
        """Average the line opacity where most of the 3x3 neighbourhood carries signal"""
        chain_outline: BooleanSocket
        """Lines between chains from a chain_id AOV pass (Composite ID Outline)"""
        chain_id: FloatSocket
        """chain_id AOV pass from the Render Layers node"""
        chain_low: FloatSocket
        """Number of neighbours in another chain, of 24, at which the line starts to appear"""
        chain_high: FloatSocket
        """Number of neighbours in another chain, of 24, at which the line is fully opaque"""
        residue_outline: BooleanSocket
        """Lines between residues from a res_id AOV pass (Composite ID Outline)"""
        residue_id: FloatSocket
        """res_id AOV pass from the Render Layers node"""
        residue_low: FloatSocket
        """Number of neighbours in another residue, of 24, at which the line starts to appear"""
        residue_high: FloatSocket
        """Number of neighbours in another residue, of 24, at which the line is fully opaque"""
        residue_difference: FloatSocket
        """Difference in residue number above which two pixels are in different residues; raise it to outline groups of residues"""

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
        base_color: InputMenu | Literal["Diffuse Color", "Image"] = "Image",
        image: InputColor = None,
        alpha: InputFloat = 1.0,
        depth: InputFloat = 0.0,
        diffuse_color: InputColor = None,
        world_scale: InputFloat = 0.1,
        shadow: InputBoolean = True,
        radius: InputFloat = 50.0,
        cone_angle: InputFloat = 2.0,
        min_gap: InputFloat = 1.0,
        strength: InputFloat = 1.0,
        max_darkening: InputFloat = 0.5,
        pixel_size: InputFloat = 0.1,
        fog: InputBoolean = False,
        near: InputFloat = 0.0,
        far: InputFloat = 100.0,
        front_fog: InputFloat = 1.0,
        back_fog: InputFloat = 1.0,
        fog_color: InputColor = None,
        contour_outline: InputBoolean = True,
        contour_high: InputFloat = 10.0,
        smooth: InputBoolean = True,
        chain_outline: InputBoolean = False,
        chain_id: InputFloat = 0.0,
        chain_low: InputFloat = 3.0,
        chain_high: InputFloat = 10.0,
        residue_outline: InputBoolean = False,
        residue_id: InputFloat = 0.0,
        residue_low: InputFloat = 3.0,
        residue_high: InputFloat = 8.0,
        residue_difference: InputFloat = 0.5,
    ):
        super().__init__(
            **{
                "Base Color": base_color,
                "Image": image,
                "Alpha": alpha,
                "Depth": depth,
                "Diffuse Color": diffuse_color,
                "World Scale": world_scale,
                "Shadow": shadow,
                "Radius": radius,
                "Cone Angle": cone_angle,
                "Min Gap": min_gap,
                "Strength": strength,
                "Max Darkening": max_darkening,
                "Pixel Size": pixel_size,
                "Fog": fog,
                "Near": near,
                "Far": far,
                "Front Fog": front_fog,
                "Back Fog": back_fog,
                "Fog Color": fog_color,
                "Contour Outline": contour_outline,
                "Contour High": contour_high,
                "Smooth": smooth,
                "Chain Outline": chain_outline,
                "Chain ID": chain_id,
                "Chain Low": chain_low,
                "Chain High": chain_high,
                "Residue Outline": residue_outline,
                "Residue ID": residue_id,
                "Residue Low": residue_low,
                "Residue High": residue_high,
                "Residue Difference": residue_difference,
            }
        )

    def _build_group(self, tree: TreeBuilder[CompositorNodeTree]) -> None:
        base_color = tree.inputs.menu(
            "Base Color",
            description="The flat Diffuse Color pass, which ignores lighting and materials, or the rendered image",
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
        diffuse_color = tree.inputs.color(
            "Diffuse Color",
            (0.8, 0.8, 0.8, 1.0),
            description="Diffuse Color pass from the Render Layers node",
            hide_value=True,
        )
        world_scale = tree.inputs.float(
            "World Scale",
            0.1,
            description="World units per Angstrom, used to convert the Angstrom inputs. Molecular Nodes imports structures at 0.1 (1 nm per world unit)",
            min_value=0.0,
            max_value=10_000.0,
        )
        with tree.inputs.panel("Shadow"):
            shadow = tree.inputs.boolean(
                "Shadow",
                True,
                description="Darken each pixel by the surrounding pixels that are closer to the camera (Composite Cone Shadow)",
                is_panel_toggle=True,
            )
            radius = tree.inputs.float(
                "Radius",
                50.0,
                description="Radius of the shadow neighbourhood in pixels",
                min_value=1.0,
                max_value=10_000.0,
            )
            cone_angle = tree.inputs.float(
                "Cone Angle",
                2.0,
                description="Tightness of the shadow cone; higher values keep shadows local to crevices",
                min_value=0.0,
                max_value=10_000.0,
            )
            min_gap = tree.inputs.float(
                "Min Gap",
                1.0,
                description="Smallest depth gap, in Angstrom, that casts a shadow",
                min_value=0.0,
                max_value=10_000.0,
            )
            strength = tree.inputs.float(
                "Strength",
                1.0,
                description="Darkening when every sample occludes the pixel",
                min_value=0.0,
                max_value=10.0,
            )
            max_darkening = tree.inputs.float(
                "Max Darkening",
                0.5,
                description="The shadow factor never drops below this value",
                min_value=0.0,
                max_value=1.0,
                subtype="FACTOR",
            )
            pixel_size = tree.inputs.float(
                "Pixel Size",
                0.1,
                description="Width of one pixel in Angstrom at the structure. canvas.compositor.illustrate() sets it from the camera",
                min_value=0.0,
                max_value=10_000.0,
            )
        with tree.inputs.panel("Fog", default_closed=True):
            fog = tree.inputs.boolean(
                "Fog",
                False,
                description="Blend towards Fog Color with depth (Composite Depth Fog)",
                is_panel_toggle=True,
            )
            near = tree.inputs.float(
                "Near",
                0.0,
                description="Depth, in world units, of the front of the structure. canvas.compositor.illustrate() sets it from the scene bounds",
                min_value=0.0,
                max_value=1_000_000.0,
            )
            far = tree.inputs.float(
                "Far",
                100.0,
                description="Depth, in world units, of the back of the structure",
                min_value=0.0,
                max_value=1_000_000.0,
            )
            front_fog = tree.inputs.float(
                "Front Fog",
                1.0,
                description="Fog factor at Near: 1 keeps the colour, 0 is fully fogged",
                min_value=0.0,
                max_value=1.0,
                subtype="FACTOR",
            )
            back_fog = tree.inputs.float(
                "Back Fog",
                1.0,
                description="Fog factor at Far: 1 keeps the colour, 0 is fully fogged",
                min_value=0.0,
                max_value=1.0,
                subtype="FACTOR",
            )
            fog_color = tree.inputs.color(
                "Fog Color",
                (1.0, 1.0, 1.0, 1.0),
                description="Colour the structure fades towards",
            )
        with tree.inputs.panel("Contour Outline"):
            contour_outline = tree.inputs.boolean(
                "Contour Outline",
                True,
                description="Lines on depth steps from a Laplacian of the depth pass (Composite Contour Outline)",
                is_panel_toggle=True,
            )
            contour_high = tree.inputs.float(
                "Contour High",
                10.0,
                description="Laplacian value, in Angstrom, at which the line is fully opaque",
                min_value=0.0,
                max_value=10_000.0,
            )
            smooth = tree.inputs.boolean(
                "Smooth",
                True,
                description="Average the line opacity where most of the 3x3 neighbourhood carries signal",
            )
        with tree.inputs.panel("Chain Outline", default_closed=True):
            chain_outline = tree.inputs.boolean(
                "Chain Outline",
                False,
                description="Lines between chains from a chain_id AOV pass (Composite ID Outline)",
                is_panel_toggle=True,
            )
            chain_id = tree.inputs.float(
                "Chain ID",
                0.0,
                description="chain_id AOV pass from the Render Layers node",
                hide_value=True,
            )
            chain_low = tree.inputs.float(
                "Chain Low",
                3.0,
                description="Number of neighbours in another chain, of 24, at which the line starts to appear",
                min_value=0.0,
                max_value=24.0,
            )
            chain_high = tree.inputs.float(
                "Chain High",
                10.0,
                description="Number of neighbours in another chain, of 24, at which the line is fully opaque",
                min_value=0.0,
                max_value=24.0,
            )
        with tree.inputs.panel("Residue Outline", default_closed=True):
            residue_outline = tree.inputs.boolean(
                "Residue Outline",
                False,
                description="Lines between residues from a res_id AOV pass (Composite ID Outline)",
                is_panel_toggle=True,
            )
            residue_id = tree.inputs.float(
                "Residue ID",
                0.0,
                description="res_id AOV pass from the Render Layers node",
                hide_value=True,
            )
            residue_low = tree.inputs.float(
                "Residue Low",
                3.0,
                description="Number of neighbours in another residue, of 24, at which the line starts to appear",
                min_value=0.0,
                max_value=24.0,
            )
            residue_high = tree.inputs.float(
                "Residue High",
                8.0,
                description="Number of neighbours in another residue, of 24, at which the line is fully opaque",
                min_value=0.0,
                max_value=24.0,
            )
            residue_difference = tree.inputs.float(
                "Residue Difference",
                0.5,
                description="Difference in residue number above which two pixels are in different residues; raise it to outline groups of residues",
                min_value=0.0,
                max_value=1_000_000.0,
            )
        image_1 = tree.outputs.color("Image", (1.0, 1.0, 1.0, 1.0))

        with c.Frame("Base colour"):
            set_alpha = c.SetAlpha(
                image=c.MenuSwitch.color(
                    base_color, {"Diffuse Color": diffuse_color, "Image": image}
                ).o.output,
                alpha=alpha,
                type="Replace Alpha",
            )
            alpha_convert = c.AlphaConvert(image=set_alpha, type="To Straight")
            _string = g.String(
                string="The base colour is given the render alpha and converted to straight alpha, so the shadow, fog and outlines act on the true colour of partly covered edge pixels. The alpha is replaced once at the end."
            )
        with c.Frame("Shadow"):
            composite_cone_shadow = CompositeConeShadow(
                depth=depth,
                radius=radius,
                cone_angle=cone_angle,
                min_gap=min_gap,
                strength=strength,
                max_darkening=max_darkening,
                pixel_size=pixel_size,
                world_scale=world_scale,
            )
            mix = g.Mix(
                a_color=alpha_convert,
                b_color=1.0 - (1.0 - composite_cone_shadow) * shadow,
                data_type="RGBA",
                blend_type="MULTIPLY",
            )
        with c.Frame("Fog"):
            composite_depth_fog = CompositeDepthFog(
                image=mix.o.result_color,
                depth=depth,
                near=near,
                far=far,
                front_fog=front_fog,
                back_fog=back_fog,
                fog_color=fog_color,
            )
            switch = fog.switch.color(mix.o.result_color, composite_depth_fog)
        with c.Frame("Outlines"):
            composite_id_outline = CompositeIDOutline(
                id=residue_id,
                min_difference=residue_difference,
                low=residue_low,
                high=residue_high,
            )
            math_1 = (
                CompositeContourOutline(
                    depth=depth, high=contour_high, smooth=smooth
                ).o.opacity
                * contour_outline
            )
            math_2 = composite_id_outline.o.opacity * residue_outline
            math_3 = (
                CompositeIDOutline(
                    id=chain_id, low=chain_low, high=chain_high
                ).o.opacity
                * chain_outline
            )
            _string_1 = g.String(
                string="Contour and residue outlines are combined by their maximum and darken the colour; the chain outline darkens it again separately. All three extend the alpha so lines survive on the background."
            )
            mix_1 = g.Mix(
                a_color=switch,
                b_color=(1.0 - math_1.max(math_2)) * (1.0 - math_3),
                data_type="RGBA",
                blend_type="MULTIPLY",
            )
            set_alpha_1 = c.SetAlpha(
                image=mix_1.o.result_color,
                alpha=alpha.max(math_1.max(math_2).max(math_3)),
                type="Replace Alpha",
            )
            alpha_convert_1 = c.AlphaConvert(image=set_alpha_1)

        alpha_convert_1 >> image_1

        base_color.default_value = "Image"


ASSET = CompositeIllustrate

ASSET_METADATA = {
    "description": "Goodsell's Illustrate composition: cone shadow, depth fog and outlines",
    "catalog_id": "441e6ca5-e514-4e77-a3cd-25fc1a2e08ae",
}
