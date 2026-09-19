# Node-group asset "Composite Depth Fog" (CompositorNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import CompositorNodeTree
from nodebpy import TreeBuilder
from nodebpy import compositor as c
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetCompositorGroup,
    ColorSocket,
    FloatSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputColor, InputFloat


class CompositeDepthFog(AssetCompositorGroup):
    """
    Linear depth fog after Goodsell's Illustrate: the image is blended towards Fog Color by a factor that runs from Front Fog at Near to Back Fog at Far (1 is no fog, 0 is fully fogged). Expects straight alpha and keeps the image alpha

    Parameters
    ----------
    image : InputColor
        Image to fog, with straight alpha
    depth : InputFloat
        Depth pass from the Render Layers node, in world units
    near : InputFloat
        Depth, in world units, of the front of the structure where the fog factor is Front Fog. canvas.compositor.illustrate() sets it from the scene bounds
    far : InputFloat
        Depth, in world units, of the back of the structure where the fog factor is Back Fog
    front_fog : InputFloat
        Fog factor at Near: 1 keeps the image colour, 0 is fully fogged
    back_fog : InputFloat
        Fog factor at Far: 1 keeps the image colour, 0 is fully fogged
    fog_color : InputColor
        Colour the image fades towards

    Inputs
    ------
    i.image : ColorSocket
        Image to fog, with straight alpha
    i.depth : FloatSocket
        Depth pass from the Render Layers node, in world units
    i.near : FloatSocket
        Depth, in world units, of the front of the structure where the fog factor is Front Fog. canvas.compositor.illustrate() sets it from the scene bounds
    i.far : FloatSocket
        Depth, in world units, of the back of the structure where the fog factor is Back Fog
    i.front_fog : FloatSocket
        Fog factor at Near: 1 keeps the image colour, 0 is fully fogged
    i.back_fog : FloatSocket
        Fog factor at Far: 1 keeps the image colour, 0 is fully fogged
    i.fog_color : ColorSocket
        Colour the image fades towards

    Outputs
    -------
    o.image : ColorSocket
        Image
    """

    _name = "Composite Depth Fog"
    _asset_name = "Composite Depth Fog"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "FILTER"
    _tree_properties = {
        "description": "Linear depth fog after Goodsell's Illustrate: the image is blended towards Fog Color by a factor that runs from Front Fog at Near to Back Fog at Far (1 is no fog, 0 is fully fogged). Expects straight alpha and keeps the image alpha"
    }

    class _Inputs(SocketAccessor):
        image: ColorSocket
        """Image to fog, with straight alpha"""
        depth: FloatSocket
        """Depth pass from the Render Layers node, in world units"""
        near: FloatSocket
        """Depth, in world units, of the front of the structure where the fog factor is Front Fog. canvas.compositor.illustrate() sets it from the scene bounds"""
        far: FloatSocket
        """Depth, in world units, of the back of the structure where the fog factor is Back Fog"""
        front_fog: FloatSocket
        """Fog factor at Near: 1 keeps the image colour, 0 is fully fogged"""
        back_fog: FloatSocket
        """Fog factor at Far: 1 keeps the image colour, 0 is fully fogged"""
        fog_color: ColorSocket
        """Colour the image fades towards"""

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
        depth: InputFloat = 0.0,
        near: InputFloat = 0.0,
        far: InputFloat = 100.0,
        front_fog: InputFloat = 1.0,
        back_fog: InputFloat = 1.0,
        fog_color: InputColor = None,
    ):
        super().__init__(
            **{
                "Image": image,
                "Depth": depth,
                "Near": near,
                "Far": far,
                "Front Fog": front_fog,
                "Back Fog": back_fog,
                "Fog Color": fog_color,
            }
        )

    def _build_group(self, tree: TreeBuilder[CompositorNodeTree]) -> None:
        image = tree.inputs.color(
            "Image",
            (0.0, 0.0, 0.0, 1.0),
            description="Image to fog, with straight alpha",
            hide_value=True,
        )
        depth = tree.inputs.float(
            "Depth",
            0.0,
            description="Depth pass from the Render Layers node, in world units",
            hide_value=True,
        )
        near = tree.inputs.float(
            "Near",
            0.0,
            description="Depth, in world units, of the front of the structure where the fog factor is Front Fog. canvas.compositor.illustrate() sets it from the scene bounds",
            min_value=0.0,
            max_value=1_000_000.0,
        )
        far = tree.inputs.float(
            "Far",
            100.0,
            description="Depth, in world units, of the back of the structure where the fog factor is Back Fog",
            min_value=0.0,
            max_value=1_000_000.0,
        )
        front_fog = tree.inputs.float(
            "Front Fog",
            1.0,
            description="Fog factor at Near: 1 keeps the image colour, 0 is fully fogged",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        back_fog = tree.inputs.float(
            "Back Fog",
            1.0,
            description="Fog factor at Far: 1 keeps the image colour, 0 is fully fogged",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        fog_color = tree.inputs.color(
            "Fog Color",
            (1.0, 1.0, 1.0, 1.0),
            description="Colour the image fades towards",
        )
        image_1 = tree.outputs.color("Image", (1.0, 1.0, 1.0, 1.0))

        with c.Frame("Fog factor"):
            mix = (
                1.0 - (front_fog - depth.map_range(near, far) * (front_fog - back_fog))
            ).mix.color(image, fog_color)
            set_alpha = c.SetAlpha(
                image=mix,
                alpha=c.SeparateColor(image=image).o.alpha,
                type="Replace Alpha",
            )
            _string = g.String(
                string="The fog factor is Front Fog at Near and Back Fog at Far, linear in depth between them and clamped outside. The colour is mixed towards Fog Color by one minus that factor, and the image alpha is kept."
            )

        set_alpha >> image_1


ASSET = CompositeDepthFog

ASSET_METADATA = {
    "description": "Linear depth fog towards a colour, after Goodsell's Illustrate",
    "catalog_id": "441e6ca5-e514-4e77-a3cd-25fc1a2e08ae",
}
