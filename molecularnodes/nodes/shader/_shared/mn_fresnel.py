# Node group "MN Fresnel" (ShaderNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import ShaderNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy import shader as s
from nodebpy.builder import CustomShaderGroup, FloatSocket, SocketAccessor
from nodebpy.types import InputFloat
from .mn_mask_transparent import MN_mask_transparent


class MNFresnel(CustomShaderGroup):
    """
    MN Fresnel

    Parameters
    ----------
    ior : InputFloat
        IOR
    factor : InputFloat
        Whether to show fresnel through transparent surfaces

    Inputs
    ------
    i.ior : FloatSocket
        IOR
    i.factor : FloatSocket
        Whether to show fresnel through transparent surfaces

    Outputs
    -------
    o.value : FloatSocket
        Value
    """

    _name = "MN Fresnel"

    class _Inputs(SocketAccessor):
        ior: FloatSocket
        """IOR"""
        factor: FloatSocket
        """Whether to show fresnel through transparent surfaces"""

    class _Outputs(SocketAccessor):
        value: FloatSocket
        """Value"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        ior: InputFloat = 0.98,
        factor: InputFloat = 0.0,
    ):
        super().__init__(**{"IOR": ior, "Factor": factor})

    def _build_group(self, tree: TreeBuilder[ShaderNodeTree]) -> None:
        ior = tree.inputs.float("IOR", 0.98, min_value=0.0, max_value=1000.0)
        factor = tree.inputs.float(
            "Factor",
            0.0,
            description="Whether to show fresnel through transparent surfaces",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        value = tree.outputs.float("Value")

        fresnel = s.Fresnel(ior=ior)
        mix = g.Mix(
            factor_float=factor,
            a_float=MN_mask_transparent(value=fresnel),
            b_float=fresnel,
            clamp_factor=True,
        )

        mix >> value
