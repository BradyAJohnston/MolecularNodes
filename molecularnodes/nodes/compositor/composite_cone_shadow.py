# Node-group asset "Composite Cone Shadow" (CompositorNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
import math
from typing import TYPE_CHECKING
from bpy.types import CompositorNodeTree
from nodebpy import TreeBuilder
from nodebpy import compositor as c
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetCompositorGroup,
    FloatSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputFloat
from ._shared.cone_shadow_sample import ConeShadowSample


class CompositeConeShadow(AssetCompositorGroup):
    """
    Conical soft shadow from the Depth pass alone, after David Goodsell's Illustrate: each pixel is darkened by the surrounding pixels that are closer to the camera by more than a gap that grows with their distance, so form comes from occlusion rather than lighting. 32 samples: 8 directions on 4 rings up to Radius

    Parameters
    ----------
    depth : InputFloat
        Depth pass from the Render Layers node, in world units
    radius : InputFloat
        Radius of the sampled neighbourhood in pixels; the four rings of samples sit at a quarter, half, three quarters and the full radius
    cone_angle : InputFloat
        Tightness of the shadow cone: a sample casts shadow only when it is closer to the camera by more than its screen distance times this slope. Higher values keep shadows local to crevices
    min_gap : InputFloat
        Smallest depth gap, in Angstrom, that casts a shadow; removes speckle in tight crevices
    strength : InputFloat
        Darkening when every sample occludes the pixel; the shadow per sample is this divided by the 32 samples
    max_darkening : InputFloat
        The shadow factor never drops below this value, so occluded pixels keep some of their colour instead of going black
    pixel_size : InputFloat
        Width of one pixel in Angstrom at the structure, which relates screen distance to depth for the cone test. canvas.compositor.illustrate() sets it from the camera
    world_scale : InputFloat
        World units per Angstrom, used to convert the Angstrom inputs. Molecular Nodes imports structures at 0.1 (1 nm per world unit)

    Inputs
    ------
    i.depth : FloatSocket
        Depth pass from the Render Layers node, in world units
    i.radius : FloatSocket
        Radius of the sampled neighbourhood in pixels; the four rings of samples sit at a quarter, half, three quarters and the full radius
    i.cone_angle : FloatSocket
        Tightness of the shadow cone: a sample casts shadow only when it is closer to the camera by more than its screen distance times this slope. Higher values keep shadows local to crevices
    i.min_gap : FloatSocket
        Smallest depth gap, in Angstrom, that casts a shadow; removes speckle in tight crevices
    i.strength : FloatSocket
        Darkening when every sample occludes the pixel; the shadow per sample is this divided by the 32 samples
    i.max_darkening : FloatSocket
        The shadow factor never drops below this value, so occluded pixels keep some of their colour instead of going black
    i.pixel_size : FloatSocket
        Width of one pixel in Angstrom at the structure, which relates screen distance to depth for the cone test. canvas.compositor.illustrate() sets it from the camera
    i.world_scale : FloatSocket
        World units per Angstrom, used to convert the Angstrom inputs. Molecular Nodes imports structures at 0.1 (1 nm per world unit)

    Outputs
    -------
    o.shadow : FloatSocket
        Shadow
    """

    _name = "Composite Cone Shadow"
    _asset_name = "Composite Cone Shadow"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "FILTER"
    _tree_properties = {
        "description": "Conical soft shadow from the Depth pass alone, after David Goodsell's Illustrate: each pixel is darkened by the surrounding pixels that are closer to the camera by more than a gap that grows with their distance, so form comes from occlusion rather than lighting. 32 samples: 8 directions on 4 rings up to Radius"
    }

    class _Inputs(SocketAccessor):
        depth: FloatSocket
        """Depth pass from the Render Layers node, in world units"""
        radius: FloatSocket
        """Radius of the sampled neighbourhood in pixels; the four rings of samples sit at a quarter, half, three quarters and the full radius"""
        cone_angle: FloatSocket
        """Tightness of the shadow cone: a sample casts shadow only when it is closer to the camera by more than its screen distance times this slope. Higher values keep shadows local to crevices"""
        min_gap: FloatSocket
        """Smallest depth gap, in Angstrom, that casts a shadow; removes speckle in tight crevices"""
        strength: FloatSocket
        """Darkening when every sample occludes the pixel; the shadow per sample is this divided by the 32 samples"""
        max_darkening: FloatSocket
        """The shadow factor never drops below this value, so occluded pixels keep some of their colour instead of going black"""
        pixel_size: FloatSocket
        """Width of one pixel in Angstrom at the structure, which relates screen distance to depth for the cone test. canvas.compositor.illustrate() sets it from the camera"""
        world_scale: FloatSocket
        """World units per Angstrom, used to convert the Angstrom inputs. Molecular Nodes imports structures at 0.1 (1 nm per world unit)"""

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
        radius: InputFloat = 50.0,
        cone_angle: InputFloat = 2.0,
        min_gap: InputFloat = 1.0,
        strength: InputFloat = 1.0,
        max_darkening: InputFloat = 0.5,
        pixel_size: InputFloat = 0.1,
        world_scale: InputFloat = 0.1,
    ):
        super().__init__(
            **{
                "Depth": depth,
                "Radius": radius,
                "Cone Angle": cone_angle,
                "Min Gap": min_gap,
                "Strength": strength,
                "Max Darkening": max_darkening,
                "Pixel Size": pixel_size,
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
        radius = tree.inputs.float(
            "Radius",
            50.0,
            description="Radius of the sampled neighbourhood in pixels; the four rings of samples sit at a quarter, half, three quarters and the full radius",
            min_value=1.0,
            max_value=10_000.0,
        )
        cone_angle = tree.inputs.float(
            "Cone Angle",
            2.0,
            description="Tightness of the shadow cone: a sample casts shadow only when it is closer to the camera by more than its screen distance times this slope. Higher values keep shadows local to crevices",
            min_value=0.0,
            max_value=10_000.0,
        )
        min_gap = tree.inputs.float(
            "Min Gap",
            1.0,
            description="Smallest depth gap, in Angstrom, that casts a shadow; removes speckle in tight crevices",
            min_value=0.0,
            max_value=10_000.0,
        )
        strength = tree.inputs.float(
            "Strength",
            1.0,
            description="Darkening when every sample occludes the pixel; the shadow per sample is this divided by the 32 samples",
            min_value=0.0,
            max_value=10.0,
        )
        max_darkening = tree.inputs.float(
            "Max Darkening",
            0.5,
            description="The shadow factor never drops below this value, so occluded pixels keep some of their colour instead of going black",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        pixel_size = tree.inputs.float(
            "Pixel Size",
            0.1,
            description="Width of one pixel in Angstrom at the structure, which relates screen distance to depth for the cone test. canvas.compositor.illustrate() sets it from the camera",
            min_value=0.0,
            max_value=10_000.0,
        )
        world_scale = tree.inputs.float(
            "World Scale",
            0.1,
            description="World units per Angstrom, used to convert the Angstrom inputs. Molecular Nodes imports structures at 0.1 (1 nm per world unit)",
            min_value=0.0,
            max_value=10_000.0,
        )
        shadow = tree.outputs.float("Shadow")

        with c.Frame("Samples"):
            math_1 = cone_angle * pixel_size * world_scale
            math_2 = min_gap * world_scale
            math_3 = radius * 0.25
            math_4 = (math_3 * math_1).max(math_2)
            math_5 = ConeShadowSample(
                depth=depth, radius=math_3, threshold=math_4
            ).o.occluded + ConeShadowSample(
                depth=depth, radius=math_3, angle=math.pi / 4, threshold=math_4
            )
            math_6 = math_5 + ConeShadowSample(
                depth=depth, radius=math_3, angle=math.pi / 2, threshold=math_4
            )
            math_7 = math_6 + ConeShadowSample(
                depth=depth, radius=math_3, angle=3 * math.pi / 4, threshold=math_4
            )
            math_8 = (
                math_7
                + ConeShadowSample(
                    depth=depth, radius=math_3, angle=math.pi, threshold=math_4
                )
                + ConeShadowSample(
                    depth=depth, radius=math_3, angle=5 * math.pi / 4, threshold=math_4
                )
            )
            math_9 = math_8 + ConeShadowSample(
                depth=depth, radius=math_3, angle=3 * math.pi / 2, threshold=math_4
            )
            math_10 = math_9 + ConeShadowSample(
                depth=depth, radius=math_3, angle=7 * math.pi / 4, threshold=math_4
            )
            math_11 = radius * 0.5
            math_12 = (math_11 * math_1).max(math_2)
            math_13 = (
                math_10
                + ConeShadowSample(depth=depth, radius=math_11, threshold=math_12)
                + ConeShadowSample(
                    depth=depth, radius=math_11, angle=math.pi / 4, threshold=math_12
                )
            )
            math_14 = math_13 + ConeShadowSample(
                depth=depth, radius=math_11, angle=math.pi / 2, threshold=math_12
            )
            math_15 = math_14 + ConeShadowSample(
                depth=depth, radius=math_11, angle=3 * math.pi / 4, threshold=math_12
            )
            math_16 = math_15 + ConeShadowSample(
                depth=depth, radius=math_11, angle=math.pi, threshold=math_12
            )
            math_17 = math_16 + ConeShadowSample(
                depth=depth, radius=math_11, angle=5 * math.pi / 4, threshold=math_12
            )
            math_18 = math_17 + ConeShadowSample(
                depth=depth, radius=math_11, angle=3 * math.pi / 2, threshold=math_12
            )
            math_19 = math_18 + ConeShadowSample(
                depth=depth, radius=math_11, angle=7 * math.pi / 4, threshold=math_12
            )
            math_20 = radius * 0.75
            math_21 = (math_20 * math_1).max(math_2)
            math_22 = (
                math_19
                + ConeShadowSample(depth=depth, radius=math_20, threshold=math_21)
                + ConeShadowSample(
                    depth=depth, radius=math_20, angle=math.pi / 4, threshold=math_21
                )
            )
            math_23 = math_22 + ConeShadowSample(
                depth=depth, radius=math_20, angle=math.pi / 2, threshold=math_21
            )
            math_24 = math_23 + ConeShadowSample(
                depth=depth, radius=math_20, angle=3 * math.pi / 4, threshold=math_21
            )
            math_25 = math_24 + ConeShadowSample(
                depth=depth, radius=math_20, angle=math.pi, threshold=math_21
            )
            math_26 = math_25 + ConeShadowSample(
                depth=depth, radius=math_20, angle=5 * math.pi / 4, threshold=math_21
            )
            math_27 = math_26 + ConeShadowSample(
                depth=depth, radius=math_20, angle=3 * math.pi / 2, threshold=math_21
            )
            math_28 = math_27 + ConeShadowSample(
                depth=depth, radius=math_20, angle=7 * math.pi / 4, threshold=math_21
            )
            math_29 = radius * 1.0
            math_30 = (math_29 * math_1).max(math_2)
            math_31 = (
                math_28
                + ConeShadowSample(depth=depth, radius=math_29, threshold=math_30)
                + ConeShadowSample(
                    depth=depth, radius=math_29, angle=math.pi / 4, threshold=math_30
                )
            )
            math_32 = math_31 + ConeShadowSample(
                depth=depth, radius=math_29, angle=math.pi / 2, threshold=math_30
            )
            math_33 = math_32 + ConeShadowSample(
                depth=depth, radius=math_29, angle=3 * math.pi / 4, threshold=math_30
            )
            math_34 = math_33 + ConeShadowSample(
                depth=depth, radius=math_29, angle=math.pi, threshold=math_30
            )
            math_35 = math_34 + ConeShadowSample(
                depth=depth, radius=math_29, angle=5 * math.pi / 4, threshold=math_30
            )
            math_36 = math_35 + ConeShadowSample(
                depth=depth, radius=math_29, angle=3 * math.pi / 2, threshold=math_30
            )
            math_37 = math_36 + ConeShadowSample(
                depth=depth, radius=math_29, angle=7 * math.pi / 4, threshold=math_30
            )
            _string = g.String(
                string="Each sample translates the depth pass by a ring radius in one of eight directions and tests whether that pixel is closer to the camera than the centre by more than the ring's threshold: the larger of Min Gap and the ring's screen distance times the cone slope, all in world units. Occluding samples are counted."
            )
        with c.Frame("Shadow"):
            (1.0 - math_37 * (strength / 32.0)).max(max_darkening) >> shadow
            _string_1 = g.String(
                string="The shadow factor starts at one and loses Strength / 32 per occluding sample, clamped at Max Darkening."
            )


ASSET = CompositeConeShadow

ASSET_METADATA = {
    "description": "Conical soft shadow from the Depth pass, after Goodsell's Illustrate",
    "catalog_id": "441e6ca5-e514-4e77-a3cd-25fc1a2e08ae",
}
