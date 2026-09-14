# Node-group asset "Transparent Outline Internal" (ShaderNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import ShaderNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy import shader as s
from nodebpy.builder import (
    AssetShaderGroup,
    ColorSocket,
    FloatSocket,
    MenuSocket,
    PackageLibrary,
    ShaderSocket,
    SocketAccessor,
)
from nodebpy.types import InputColor, InputFloat, InputMenu
from ._shared.mn_fresnel import MNFresnel
from .mn_color import MNColor


class TransparentOutlineInternal(AssetShaderGroup):
    """
    Transparent Outline Internal

    Parameters
    ----------
    transparency : InputFloat
        Blend weight to use for mixing two shaders. At zero it uses the first shader entirely and at one the second shader
    menu : InputMenu | Literal["Transparent", "Fresnel"]
        Menu
    outline_color : InputColor
        Outline Color

    Inputs
    ------
    i.transparency : FloatSocket
        Blend weight to use for mixing two shaders. At zero it uses the first shader entirely and at one the second shader
    i.menu : MenuSocket
        Menu
    i.outline_color : ColorSocket
        Outline Color

    Outputs
    -------
    o.shader : ShaderSocket
        Shader
    """

    _name = "Transparent Outline Internal"
    _asset_name = "Transparent Outline Internal"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")

    class _Inputs(SocketAccessor):
        transparency: FloatSocket
        """Blend weight to use for mixing two shaders. At zero it uses the first shader entirely and at one the second shader"""
        menu: MenuSocket
        """Menu"""
        outline_color: ColorSocket
        """Outline Color"""

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
        transparency: InputFloat = 0.9,
        menu: InputMenu | Literal["Transparent", "Fresnel"] = "Transparent",
        outline_color: InputColor = None,
    ):
        super().__init__(
            **{
                "Transparency": transparency,
                "Menu": menu,
                "Outline Color": outline_color,
            }
        )

    def _build_group(self, tree: TreeBuilder[ShaderNodeTree]) -> None:
        transparency = tree.inputs.float(
            "Transparency",
            0.9,
            description="Blend weight to use for mixing two shaders. At zero it uses the first shader entirely and at one the second shader",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        menu = tree.inputs.menu("Menu", expanded=True, optional_label=True)
        outline_color = tree.inputs.color("Outline Color", (1.0, 1.0, 1.0, 1.0))
        shader = tree.outputs.shader("Shader")

        mix_shader = s.MixShader(
            fac=g.Math.greater_than(s.LightPath().o.transparent_depth, 0.0).o.value
            + transparency,
            shader=s.DiffuseBSDF(color=MNColor().o.color),
            shader_001=s.TransparentBSDF(),
        )
        (
            s.MenuSwitch.shader(
                menu,
                {
                    "Transparent": mix_shader,
                    "Fresnel": s.MixShader(
                        fac=MNFresnel(ior=0.95),
                        shader=mix_shader,
                        shader_001=outline_color,
                    ),
                },
            )
            >> shader
        )

        menu.default_value = "Transparent"


ASSET = TransparentOutlineInternal

ASSET_METADATA = {
    "catalog_id": "fc8d3698-34f7-4b7e-8167-a2c0391b171b",
}
