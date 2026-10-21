# Node-group asset "Cone Shadow" (CompositorNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
import math
from typing import TYPE_CHECKING
from bpy.types import CompositorNodeTree
from nodebpy import TreeBuilder
from nodebpy import compositor as c
from nodebpy.builder import (
    AssetCompositorGroup,
    CustomCompositorGroup,
    FloatSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputFloat


class ConeShadowSample(CustomCompositorGroup):
    _name = "Cone Shadow Sample"
    _color_tag = "FILTER"
    _tree_properties = {
        "description": "One sample of the conical shadow: whether the depth at a pixel offset is closer to the camera than the pixel by more than a threshold"
    }

    def _build_group(self, tree: TreeBuilder[CompositorNodeTree]) -> None:
        depth = tree.inputs.float(
            "Depth", 0.0, description="Depth pass, in world units", hide_value=True
        )
        radius = tree.inputs.float(
            "Radius",
            0.0,
            description="Distance of the sample from the pixel, in pixels",
            min_value=0.0,
            max_value=10_000.0,
        )
        angle = tree.inputs.float(
            "Angle",
            0.0,
            description="Direction of the sample from the pixel",
            subtype="ANGLE",
        )
        threshold = tree.inputs.float(
            "Threshold",
            0.0,
            description="Depth gap, in world units, above which the sample casts shadow on the pixel",
            min_value=0.0,
            max_value=10_000.0,
        )
        occluded = tree.outputs.float("Occluded")

        math_1 = depth.min(10_000.0)
        translate = c.Translate(
            image=math_1,
            x=radius * angle.cos(),
            y=radius * angle.sin(),
            interpolation="Nearest",
            extension_x="Extend",
            extension_y="Extend",
        )
        math_2 = c.Math.greater_than(math_1 - translate, threshold)

        math_2 >> occluded


class ConeShadow(AssetCompositorGroup):
    """
    Conical soft shadow from the Depth pass alone, after David Goodsell's Illustrate: each pixel is darkened by the surrounding pixels that are closer to the camera by more than a gap that grows with their distance, so form comes from occlusion rather than lighting. 32 samples: 8 directions on 4 rings up to Radius

    Parameters
    ----------
    depth : InputFloat
        Depth pass from the Render Layers node, in world units
    pixel_radius : InputFloat
        Radius of the sampled neighbourhood in pixels; the four rings of samples sit at a quarter, half, three quarters and the full radius
    pixel_size : InputFloat
        Width of one pixel in world units at the structure, like the depth pass, which relates screen distance to depth for the cone test. canvas.compositor.illustrate() sets it from the camera
    angle : InputFloat
        Tightness of the shadow cone: a sample casts shadow only when it is closer to the camera by more than its screen distance times this slope. Higher values keep shadows local to crevices
    threshold : InputFloat
        Smallest depth gap, in Nanometres, that casts a shadow; removes speckle in tight crevices
    shadow_strength : InputFloat
        Darkening when every sample occludes the pixel; the shadow per sample is this divided by the 32 samples
    shadow_max : InputFloat
        The shadow factor never drops below this value, so occluded pixels keep some of their colour instead of going black

    Inputs
    ------
    i.depth : FloatSocket
        Depth pass from the Render Layers node, in world units
    i.pixel_radius : FloatSocket
        Radius of the sampled neighbourhood in pixels; the four rings of samples sit at a quarter, half, three quarters and the full radius
    i.pixel_size : FloatSocket
        Width of one pixel in world units at the structure, like the depth pass, which relates screen distance to depth for the cone test. canvas.compositor.illustrate() sets it from the camera
    i.angle : FloatSocket
        Tightness of the shadow cone: a sample casts shadow only when it is closer to the camera by more than its screen distance times this slope. Higher values keep shadows local to crevices
    i.threshold : FloatSocket
        Smallest depth gap, in Nanometres, that casts a shadow; removes speckle in tight crevices
    i.shadow_strength : FloatSocket
        Darkening when every sample occludes the pixel; the shadow per sample is this divided by the 32 samples
    i.shadow_max : FloatSocket
        The shadow factor never drops below this value, so occluded pixels keep some of their colour instead of going black

    Outputs
    -------
    o.shadow : FloatSocket
        Shadow
    """

    _name = "Cone Shadow"
    _asset_name = "Cone Shadow"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "FILTER"
    _tree_properties = {
        "description": "Conical soft shadow from the Depth pass alone, after David Goodsell's Illustrate: each pixel is darkened by the surrounding pixels that are closer to the camera by more than a gap that grows with their distance, so form comes from occlusion rather than lighting. 32 samples: 8 directions on 4 rings up to Radius"
    }

    class _Inputs(SocketAccessor):
        depth: FloatSocket
        """Depth pass from the Render Layers node, in world units"""
        pixel_radius: FloatSocket
        """Radius of the sampled neighbourhood in pixels; the four rings of samples sit at a quarter, half, three quarters and the full radius"""
        pixel_size: FloatSocket
        """Width of one pixel in world units at the structure, like the depth pass, which relates screen distance to depth for the cone test. canvas.compositor.illustrate() sets it from the camera"""
        angle: FloatSocket
        """Tightness of the shadow cone: a sample casts shadow only when it is closer to the camera by more than its screen distance times this slope. Higher values keep shadows local to crevices"""
        threshold: FloatSocket
        """Smallest depth gap, in Nanometres, that casts a shadow; removes speckle in tight crevices"""
        shadow_strength: FloatSocket
        """Darkening when every sample occludes the pixel; the shadow per sample is this divided by the 32 samples"""
        shadow_max: FloatSocket
        """The shadow factor never drops below this value, so occluded pixels keep some of their colour instead of going black"""

    class _Outputs(SocketAccessor):
        shadow: FloatSocket
        """Shadow"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        depth: InputFloat = 0.0,
        pixel_radius: InputFloat = 50.0,
        pixel_size: InputFloat = 0.1,
        angle: InputFloat = 0.4363323,
        threshold: InputFloat = 2.0,
        shadow_strength: InputFloat = 1.0,
        shadow_max: InputFloat = 0.3,
    ):
        super().__init__(
            **{
                "Depth": depth,
                "Pixel Radius": pixel_radius,
                "Pixel Size": pixel_size,
                "Angle": angle,
                "Threshold": threshold,
                "Shadow Strength": shadow_strength,
                "Shadow Max": shadow_max,
            }
        )

    def _build_group(self, tree: TreeBuilder[CompositorNodeTree]) -> None:
        depth = tree.inputs.float(
            "Depth",
            0.0,
            description="Depth pass from the Render Layers node, in world units",
            hide_value=True,
        )
        with tree.inputs.panel("Pixel"):
            pixel_radius = tree.inputs.float(
                "Pixel Radius",
                50.0,
                description="Radius of the sampled neighbourhood in pixels; the four rings of samples sit at a quarter, half, three quarters and the full radius",
                min_value=1.0,
                max_value=10_000.0,
                subtype="PIXEL",
            )
            pixel_size = tree.inputs.float(
                "Pixel Size",
                0.1,
                description="Width of one pixel in world units at the structure, like the depth pass, which relates screen distance to depth for the cone test. canvas.compositor.illustrate() sets it from the camera",
                min_value=0.0,
                max_value=10_000.0,
            )
        with tree.inputs.panel("Cone"):
            angle = tree.inputs.float(
                "Angle",
                0.4363323,
                description="Tightness of the shadow cone: a sample casts shadow only when it is closer to the camera by more than its screen distance times this slope. Higher values keep shadows local to crevices",
                min_value=0.0,
                max_value=math.pi / 2,
                subtype="ANGLE",
            )
            threshold = tree.inputs.float(
                "Threshold",
                2.0,
                description="Smallest depth gap, in Nanometres, that casts a shadow; removes speckle in tight crevices",
                min_value=0.0,
                max_value=10_000.0,
                subtype="DISTANCE",
            )
        with tree.inputs.panel("Shadow"):
            shadow_strength = tree.inputs.float(
                "Shadow Strength",
                1.0,
                description="Darkening when every sample occludes the pixel; the shadow per sample is this divided by the 32 samples",
                min_value=0.0,
                max_value=1.0,
                subtype="FACTOR",
            )
            shadow_max = tree.inputs.float(
                "Shadow Max",
                0.3,
                description="The shadow factor never drops below this value, so occluded pixels keep some of their colour instead of going black",
                min_value=0.0,
                max_value=1.0,
                subtype="FACTOR",
            )
        shadow = tree.outputs.float("Shadow")

        with c.Frame("Samples"):
            math_1 = angle * pixel_size
            math_2 = pixel_radius * 0.25
            math_3 = (math_2 * math_1).max(threshold)
            math_4 = ConeShadowSample(
                Depth=depth, Radius=math_2, Threshold=math_3
            ).o.occluded + ConeShadowSample(
                Depth=depth, Radius=math_2, Angle=math.pi / 4, Threshold=math_3
            )
            math_5 = math_4 + ConeShadowSample(
                Depth=depth, Radius=math_2, Angle=math.pi / 2, Threshold=math_3
            )
            math_6 = math_5 + ConeShadowSample(
                Depth=depth, Radius=math_2, Angle=3 * math.pi / 4, Threshold=math_3
            )
            math_7 = (
                math_6
                + ConeShadowSample(
                    Depth=depth, Radius=math_2, Angle=math.pi, Threshold=math_3
                )
                + ConeShadowSample(
                    Depth=depth, Radius=math_2, Angle=5 * math.pi / 4, Threshold=math_3
                )
            )
            math_8 = math_7 + ConeShadowSample(
                Depth=depth, Radius=math_2, Angle=3 * math.pi / 2, Threshold=math_3
            )
            math_9 = math_8 + ConeShadowSample(
                Depth=depth, Radius=math_2, Angle=7 * math.pi / 4, Threshold=math_3
            )
            math_10 = pixel_radius * 0.5
            math_11 = (math_10 * math_1).max(threshold)
            math_12 = (
                math_9
                + ConeShadowSample(Depth=depth, Radius=math_10, Threshold=math_11)
                + ConeShadowSample(
                    Depth=depth, Radius=math_10, Angle=math.pi / 4, Threshold=math_11
                )
            )
            math_13 = math_12 + ConeShadowSample(
                Depth=depth, Radius=math_10, Angle=math.pi / 2, Threshold=math_11
            )
            math_14 = math_13 + ConeShadowSample(
                Depth=depth, Radius=math_10, Angle=3 * math.pi / 4, Threshold=math_11
            )
            math_15 = math_14 + ConeShadowSample(
                Depth=depth, Radius=math_10, Angle=math.pi, Threshold=math_11
            )
            math_16 = math_15 + ConeShadowSample(
                Depth=depth, Radius=math_10, Angle=5 * math.pi / 4, Threshold=math_11
            )
            math_17 = math_16 + ConeShadowSample(
                Depth=depth, Radius=math_10, Angle=3 * math.pi / 2, Threshold=math_11
            )
            math_18 = math_17 + ConeShadowSample(
                Depth=depth, Radius=math_10, Angle=7 * math.pi / 4, Threshold=math_11
            )
            math_19 = pixel_radius * 0.75
            math_20 = (math_19 * math_1).max(threshold)
            math_21 = (
                math_18
                + ConeShadowSample(Depth=depth, Radius=math_19, Threshold=math_20)
                + ConeShadowSample(
                    Depth=depth, Radius=math_19, Angle=math.pi / 4, Threshold=math_20
                )
            )
            math_22 = math_21 + ConeShadowSample(
                Depth=depth, Radius=math_19, Angle=math.pi / 2, Threshold=math_20
            )
            math_23 = math_22 + ConeShadowSample(
                Depth=depth, Radius=math_19, Angle=3 * math.pi / 4, Threshold=math_20
            )
            math_24 = math_23 + ConeShadowSample(
                Depth=depth, Radius=math_19, Angle=math.pi, Threshold=math_20
            )
            math_25 = math_24 + ConeShadowSample(
                Depth=depth, Radius=math_19, Angle=5 * math.pi / 4, Threshold=math_20
            )
            math_26 = math_25 + ConeShadowSample(
                Depth=depth, Radius=math_19, Angle=3 * math.pi / 2, Threshold=math_20
            )
            math_27 = math_26 + ConeShadowSample(
                Depth=depth, Radius=math_19, Angle=7 * math.pi / 4, Threshold=math_20
            )
            math_28 = pixel_radius * 1.0
            math_29 = (math_28 * math_1).max(threshold)
            math_30 = (
                math_27
                + ConeShadowSample(Depth=depth, Radius=math_28, Threshold=math_29)
                + ConeShadowSample(
                    Depth=depth, Radius=math_28, Angle=math.pi / 4, Threshold=math_29
                )
            )
            math_31 = math_30 + ConeShadowSample(
                Depth=depth, Radius=math_28, Angle=math.pi / 2, Threshold=math_29
            )
            math_32 = math_31 + ConeShadowSample(
                Depth=depth, Radius=math_28, Angle=3 * math.pi / 4, Threshold=math_29
            )
            math_33 = math_32 + ConeShadowSample(
                Depth=depth, Radius=math_28, Angle=math.pi, Threshold=math_29
            )
            math_34 = math_33 + ConeShadowSample(
                Depth=depth, Radius=math_28, Angle=5 * math.pi / 4, Threshold=math_29
            )
            math_35 = math_34 + ConeShadowSample(
                Depth=depth, Radius=math_28, Angle=3 * math.pi / 2, Threshold=math_29
            )
            math_36 = math_35 + ConeShadowSample(
                Depth=depth, Radius=math_28, Angle=7 * math.pi / 4, Threshold=math_29
            )
            _string = c.String(
                "Each sample translates the depth pass by a ring radius in one of eight directions and tests whether that pixel is closer to the camera than the centre by more than the ring's threshold: the larger of Min Gap and the ring's screen distance times the cone slope, all in world units. Occluding samples are counted."
            )
        with c.Frame("Shadow"):
            (1.0 - math_36 * (shadow_strength / 32.0)).max(shadow_max) >> shadow
            _string_1 = c.String(
                "The shadow factor starts at one and loses Strength / 32 per occluding sample, clamped at Max Darkening."
            )


ASSET = ConeShadow

ASSET_METADATA = {
    "description": "Conical soft shadow from the Depth pass, after Goodsell's Illustrate",
    "catalog_id": "874a72a6-7cf2-456a-ab07-e79f5a1e4294",
}
