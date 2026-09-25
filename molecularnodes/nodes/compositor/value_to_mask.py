# Node-group asset "Value to Mask" (CompositorNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import CompositorNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetCompositorGroup,
    FloatSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputFloat


class ValueToMask(AssetCompositorGroup):
    """
    Mask of the pixels where a float pass, such as an AOV written by a material, matches a value within a tolerance

    Parameters
    ----------
    value : InputFloat
        A float pass, such as an AOV written by a material from an attribute like chain_id
    target : InputFloat
        The value to select from the pass
    epsilon : InputFloat
        Tolerance around Target that still counts as a match

    Inputs
    ------
    i.value : FloatSocket
        A float pass, such as an AOV written by a material from an attribute like chain_id
    i.target : FloatSocket
        The value to select from the pass
    i.epsilon : FloatSocket
        Tolerance around Target that still counts as a match

    Outputs
    -------
    o.mask : FloatSocket
        Mask
    o.inverse : FloatSocket
        Inverse
    """

    _name = "Value to Mask"
    _asset_name = "Value to Mask"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {
        "description": "Mask of the pixels where a float pass, such as an AOV written by a material, matches a value within a tolerance"
    }

    class _Inputs(SocketAccessor):
        value: FloatSocket
        """A float pass, such as an AOV written by a material from an attribute like chain_id"""
        target: FloatSocket
        """The value to select from the pass"""
        epsilon: FloatSocket
        """Tolerance around Target that still counts as a match"""

    class _Outputs(SocketAccessor):
        mask: FloatSocket
        """Mask"""
        inverse: FloatSocket
        """Inverse"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        value: InputFloat = 0.0,
        target: InputFloat = 0.0,
        epsilon: InputFloat = 0.01,
    ):
        super().__init__(**{"Value": value, "Target": target, "Epsilon": epsilon})

    def _build_group(self, tree: TreeBuilder[CompositorNodeTree]) -> None:
        value = tree.inputs.float(
            "Value",
            0.0,
            description="A float pass, such as an AOV written by a material from an attribute like chain_id",
            hide_value=True,
        )
        target = tree.inputs.float(
            "Target", 0.0, description="The value to select from the pass"
        )
        epsilon = tree.inputs.float(
            "Epsilon",
            0.01,
            description="Tolerance around Target that still counts as a match",
            min_value=0.0,
            max_value=10_000.0,
        )
        mask = tree.outputs.float("Mask")
        inverse = tree.outputs.float("Inverse")

        math_1 = g.Math.compare(value, target, epsilon)
        1.0 - math_1 >> inverse

        math_1 >> mask


ASSET = ValueToMask

ASSET_METADATA = {
    "description": "Mask of the pixels where a float pass matches a value",
    "catalog_id": "441e6ca5-e514-4e77-a3cd-25fc1a2e08ae",
}
