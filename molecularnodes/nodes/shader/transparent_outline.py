# Node-group asset 'Transparent Outline' (ShaderNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from nodebpy import geometry as g
from nodebpy import shader as s
from nodebpy.builder import (
    AssetShaderGroup,
    ColorSocket,
    CustomShaderGroup,
    FloatSocket,
    MenuSocket,
    PackageLibrary,
    ShaderSocket,
    SocketAccessor,
)
from nodebpy.types import InputColor, InputFloat, InputMenu
from .mn_color import MNColor
from .outline_mask import OutlineMask


class MN_mask_transparent(CustomShaderGroup):
    _name = ".MN_mask_transparent"

    def _build_group(self, tree):
        value = tree.inputs.float("Value", 0.0, min_value=-10_000.0, max_value=10_000.0)
        value_1 = tree.outputs.float("Value")

        (
            g.Math.less_than(s.LightPath().o.transparent_depth, 1.0).o.value * value
            >> value_1
        )


class MNFresnel(CustomShaderGroup):
    _name = "MN Fresnel"

    def _build_group(self, tree):
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
            a_float=MN_mask_transparent(Value=fresnel),
            b_float=fresnel,
            clamp_factor=True,
        )

        mix >> value


class TransparentOutline(AssetShaderGroup):
    """
    Transparent Outline

    Parameters
    ----------
    alpha : InputFloat
        Blend weight to use for mixing two shaders. At zero it uses the first shader entirely and at one the second shader
    menu : InputMenu | Literal["Transparent", "Outline"]
        Menu
    outline_color : InputColor
        Outline Color
    threshold : InputFloat
        Threshold
    thickness : InputFloat
        Thickness

    Inputs
    ------
    i.alpha : FloatSocket
        Blend weight to use for mixing two shaders. At zero it uses the first shader entirely and at one the second shader
    i.menu : MenuSocket
        Menu
    i.outline_color : ColorSocket
        Outline Color
    i.threshold : FloatSocket
        Threshold
    i.thickness : FloatSocket
        Thickness

    Outputs
    -------
    o.shader : ShaderSocket
        Shader
    """

    _name = "Transparent Outline"
    _asset_name = "Transparent Outline"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")

    class _Inputs(SocketAccessor):
        alpha: FloatSocket
        """Blend weight to use for mixing two shaders. At zero it uses the first shader entirely and at one the second shader"""
        menu: MenuSocket
        """Menu"""
        outline_color: ColorSocket
        """Outline Color"""
        threshold: FloatSocket
        """Threshold"""
        thickness: FloatSocket
        """Thickness"""

    class _Outputs(SocketAccessor):
        shader: ShaderSocket
        """Shader"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        alpha: InputFloat = 0.95,
        menu: InputMenu | Literal["Transparent", "Outline"] = "Transparent",
        outline_color: InputColor = None,
        threshold: InputFloat = 0.2,
        thickness: InputFloat = 0.15,
    ):
        super().__init__(
            **{
                "Alpha": alpha,
                "Menu": menu,
                "Outline Color": outline_color,
                "Threshold": threshold,
                "Thickness": thickness,
            }
        )

    def _build_group(self, tree):
        alpha = tree.inputs.float(
            "Alpha",
            0.95,
            description="Blend weight to use for mixing two shaders. At zero it uses the first shader entirely and at one the second shader",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        menu = tree.inputs.menu("Menu", expanded=True, optional_label=True)
        outline_color = tree.inputs.color("Outline Color", (1.0, 1.0, 1.0, 1.0))
        threshold = tree.inputs.float(
            "Threshold", 0.2, min_value=0.0, max_value=10_000.0
        )
        thickness = tree.inputs.float(
            "Thickness", 0.15, min_value=0.0, max_value=10_000.0
        )
        shader = tree.outputs.shader("Shader")

        _group = MNFresnel(IOR=0.95)
        mix_shader = s.MixShader(
            fac=g.Math.greater_than(s.LightPath().o.transparent_depth, 0.0).o.value
            + alpha,
            shader=s.DiffuseBSDF(color=MNColor().o.color),
            shader_001=s.TransparentBSDF(),
        )
        mix_shader_1 = s.MixShader(
            fac=OutlineMask(threshold=threshold, thickness=thickness),
            shader=mix_shader,
            shader_001=outline_color,
        )
        (
            s.MenuSwitch.shader(
                menu, {"Transparent": mix_shader, "Outline": mix_shader_1}
            )
            >> shader
        )

        menu.default_value = "Transparent"


ASSET = TransparentOutline

ASSET_METADATA = {
    "catalog_id": "fc8d3698-34f7-4b7e-8167-a2c0391b171b",
}
