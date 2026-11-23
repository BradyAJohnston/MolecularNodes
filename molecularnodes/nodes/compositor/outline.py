# Node-group asset "Outline" (CompositorNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import CompositorNodeTree
from nodebpy import TreeBuilder
from nodebpy import compositor as c
from nodebpy.builder import (
    AssetCompositorGroup,
    ColorSocket,
    FloatSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputColor, InputFloat


class Outline(AssetCompositorGroup):
    """
    Draw a line mask over an image in a flat colour, keeping the lines visible over a transparent background

    Parameters
    ----------
    image : InputColor
        Image to draw the lines over
    mask : InputFloat
        Where to draw lines, such as the output of Composite Outline Mask
    color : InputColor
        Colour of the lines

    Inputs
    ------
    i.image : ColorSocket
        Image to draw the lines over
    i.mask : FloatSocket
        Where to draw lines, such as the output of Composite Outline Mask
    i.color : ColorSocket
        Colour of the lines

    Outputs
    -------
    o.image : ColorSocket
        Image
    """

    _name = "Outline"
    _asset_name = "Outline"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "FILTER"
    _tree_properties = {
        "description": "Draw a line mask over an image in a flat colour, keeping the lines visible over a transparent background"
    }

    class _Inputs(SocketAccessor):
        image: ColorSocket
        """Image to draw the lines over"""
        mask: FloatSocket
        """Where to draw lines, such as the output of Composite Outline Mask"""
        color: ColorSocket
        """Colour of the lines"""

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
        mask: InputFloat = 0.0,
        color: InputColor = None,
    ):
        super().__init__(**{"Image": image, "Mask": mask, "Color": color})

    def _build_group(self, tree: TreeBuilder[CompositorNodeTree]) -> None:
        image = tree.inputs.color(
            "Image",
            (0.0, 0.0, 0.0, 1.0),
            description="Image to draw the lines over",
            hide_value=True,
        )
        mask = tree.inputs.float(
            "Mask",
            0.0,
            description="Where to draw lines, such as the output of Composite Outline Mask",
            min_value=0.0,
            max_value=1.0,
            hide_value=True,
            subtype="FACTOR",
        )
        color = tree.inputs.color(
            "Color", (0.0, 0.0, 0.0, 1.0), description="Colour of the lines"
        )
        image_1 = tree.outputs.color("Image", (1.0, 1.0, 1.0, 1.0))

        with c.Frame("Draw lines"):
            alpha_over = c.AlphaOver(background=image, foreground=color, factor=mask)
        with c.Frame("Restore alpha"):
            set_alpha = c.SetAlpha(
                image=alpha_over,
                alpha=(mask * c.SeparateColor(color).o.alpha).max(
                    c.SeparateColor(image).o.alpha
                ),
                type="Replace Alpha",
            )
            _string = c.String(
                "Lines drawn over a transparent background would otherwise stay transparent. The alpha is set once, deliberately, to the larger of the line coverage and the image alpha, so lines survive on transparent backgrounds and the image alpha is kept everywhere else."
            )

        set_alpha >> image_1


ASSET = Outline

ASSET_METADATA = {
    "description": "Draw a line mask over an image in a flat colour",
    "catalog_id": "874a72a6-7cf2-456a-ab07-e79f5a1e4294",
}
