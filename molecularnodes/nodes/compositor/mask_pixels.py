# Node-group asset "Mask Pixels" (CompositorNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
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


class MaskPixels(AssetCompositorGroup):
    """
    Split an image by a mask into the pixels inside it and the pixels outside it, each with the rest made transparent, so the two can be processed separately and recombined with Alpha Over

    Parameters
    ----------
    image : InputColor
        Image to split
    mask : InputFloat
        Pixels to keep in the Masked output, such as from Value to Mask

    Inputs
    ------
    i.image : ColorSocket
        Image to split
    i.mask : FloatSocket
        Pixels to keep in the Masked output, such as from Value to Mask

    Outputs
    -------
    o.masked : ColorSocket
        Masked
    o.inverse : ColorSocket
        Inverse
    """

    _name = "Mask Pixels"
    _asset_name = "Mask Pixels"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {
        "description": "Split an image by a mask into the pixels inside it and the pixels outside it, each with the rest made transparent, so the two can be processed separately and recombined with Alpha Over"
    }

    class _Inputs(SocketAccessor):
        image: ColorSocket
        """Image to split"""
        mask: FloatSocket
        """Pixels to keep in the Masked output, such as from Value to Mask"""

    class _Outputs(SocketAccessor):
        masked: ColorSocket
        """Masked"""
        inverse: ColorSocket
        """Inverse"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        image: InputColor = None,
        mask: InputFloat = 1.0,
    ):
        super().__init__(**{"Image": image, "Mask": mask})

    def _build_group(self, tree: TreeBuilder[CompositorNodeTree]) -> None:
        image = tree.inputs.color(
            "Image", (0.0, 0.0, 0.0, 1.0), description="Image to split", hide_value=True
        )
        mask = tree.inputs.float(
            "Mask",
            1.0,
            description="Pixels to keep in the Masked output, such as from Value to Mask",
            min_value=0.0,
            max_value=1.0,
            hide_value=True,
            subtype="FACTOR",
        )
        masked = tree.outputs.color("Masked", (1.0, 1.0, 1.0, 1.0))
        inverse = tree.outputs.color("Inverse", (1.0, 1.0, 1.0, 1.0))

        set_alpha = c.SetAlpha(image=image, alpha=mask)
        set_alpha_1 = c.SetAlpha(image=image, alpha=1.0 - mask)

        set_alpha >> masked
        set_alpha_1 >> inverse


ASSET = MaskPixels

ASSET_METADATA = {
    "description": "Split an image by a mask into the pixels inside and outside it",
    "catalog_id": "441e6ca5-e514-4e77-a3cd-25fc1a2e08ae",
}
