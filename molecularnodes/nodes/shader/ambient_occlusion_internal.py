# Node-group asset "Ambient Occlusion Internal" (ShaderNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import ShaderNodeTree
from nodebpy import TreeBuilder
from nodebpy import shader as s
from nodebpy.builder import (
    AssetShaderGroup,
    FloatSocket,
    MenuSocket,
    PackageLibrary,
    ShaderSocket,
    SocketAccessor,
)
from nodebpy.types import InputFloat, InputMenu
from .color_ao import ColorAO
from .mn_color import MNColor


class AmbientOcclusionInternal(AssetShaderGroup):
    """
    Ambient Occlusion Internal

    Parameters
    ----------
    menu : InputMenu | Literal["AO", "None"]
        Shade with ambient occlusion (`AO`) or with the plain color (`None`)
    ao_space : InputMenu | Literal["Global", "Local"]
        Look in local geometry or world space for AO calculations
    distance : InputFloat
        Distance for AO calculations
    exponent : InputFloat
        Exponent to apply to AO calculations

    Inputs
    ------
    i.menu : MenuSocket
        Shade with ambient occlusion (`AO`) or with the plain color (`None`)
    i.ao_space : MenuSocket
        Look in local geometry or world space for AO calculations
    i.distance : FloatSocket
        Distance for AO calculations
    i.exponent : FloatSocket
        Exponent to apply to AO calculations

    Outputs
    -------
    o.shader : ShaderSocket
        Emission of the occluded color, mixed towards transparent by the alpha of the `Color` attribute
    """

    _name = "Ambient Occlusion Internal"
    _asset_name = "Ambient Occlusion Internal"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "SHADER"

    class _Inputs(SocketAccessor):
        menu: MenuSocket
        """Shade with ambient occlusion (`AO`) or with the plain color (`None`)"""
        ao_space: MenuSocket
        """Look in local geometry or world space for AO calculations"""
        distance: FloatSocket
        """Distance for AO calculations"""
        exponent: FloatSocket
        """Exponent to apply to AO calculations"""

    class _Outputs(SocketAccessor):
        shader: ShaderSocket
        """Emission of the occluded color, mixed towards transparent by the alpha of the `Color` attribute"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        menu: InputMenu | Literal["AO", "None"] = "AO",
        ao_space: InputMenu | Literal["Global", "Local"] = "Global",
        distance: InputFloat = 1.0,
        exponent: InputFloat = 2.0,
    ):
        super().__init__(
            **{
                "Menu": menu,
                "AO Space": ao_space,
                "Distance": distance,
                "Exponent": exponent,
            }
        )

    def _build_group(self, tree: TreeBuilder[ShaderNodeTree]) -> None:
        menu = tree.inputs.menu(
            "Menu",
            description="Shade with ambient occlusion (`AO`) or with the plain color (`None`)",
            expanded=True,
            optional_label=True,
        )
        ao_space = tree.inputs.menu(
            "AO Space",
            description="Look in local geometry or world space for AO calculations",
            expanded=True,
            optional_label=True,
        )
        distance = tree.inputs.float(
            "Distance",
            1.0,
            description="Distance for AO calculations",
            min_value=0.0,
            max_value=1000.0,
        )
        exponent = tree.inputs.float(
            "Exponent",
            2.0,
            description="Exponent to apply to AO calculations",
            min_value=0.0,
            max_value=10_000.0,
        )
        shader = tree.outputs.shader(
            "Shader",
            description="Emission of the occluded color, mixed towards transparent by the alpha of the `Color` attribute",
        )

        mn_color = MNColor()
        color_ao = ColorAO(
            color=mn_color.o.color,
            menu=menu,
            ao_space=ao_space,
            distance=distance,
            exponent=exponent,
        )
        mix_shader = s.MixShader(
            fac=mn_color.o.alpha,
            shader=s.TransparentBSDF(),
            shader_001=s.Emission(color=color_ao),
        )

        mix_shader >> shader

        menu.default_value = "AO"
        ao_space.default_value = "Global"


ASSET = AmbientOcclusionInternal

ASSET_METADATA = {
    "catalog_id": "fc8d3698-34f7-4b7e-8167-a2c0391b171b",
}
