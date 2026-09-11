# Node-group asset "Flat" (ShaderNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import ShaderNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
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
from .mn_color import MNColor
from .outline_mask import OutlineMask


class Flat(AssetShaderGroup):
    """
    Flat

    Parameters
    ----------
    outline : InputMenu | Literal["Outline", "None"]
        Outline
    threshold : InputFloat
        Threshold
    thickness : InputFloat
        Thickness

    Inputs
    ------
    i.outline : MenuSocket
        Outline
    i.threshold : FloatSocket
        Threshold
    i.thickness : FloatSocket
        Thickness

    Outputs
    -------
    o.emission : ShaderSocket
        Emission
    """

    _name = "Flat"
    _asset_name = "Flat"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "SHADER"

    class _Inputs(SocketAccessor):
        outline: MenuSocket
        """Outline"""
        threshold: FloatSocket
        """Threshold"""
        thickness: FloatSocket
        """Thickness"""

    class _Outputs(SocketAccessor):
        emission: ShaderSocket
        """Emission"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        outline: InputMenu | Literal["Outline", "None"] = "Outline",
        threshold: InputFloat = 0.8,
        thickness: InputFloat = 0.15,
    ):
        super().__init__(
            **{"Outline": outline, "Threshold": threshold, "Thickness": thickness}
        )

    def _build_group(self, tree: TreeBuilder[ShaderNodeTree]) -> None:
        outline = tree.inputs.menu("Outline", expanded=True, optional_label=True)
        threshold = tree.inputs.float(
            "Threshold", 0.8, min_value=0.0, max_value=10_000.0
        )
        thickness = tree.inputs.float(
            "Thickness", 0.15, min_value=0.0, max_value=10_000.0
        )
        emission = tree.outputs.shader("Emission")

        group = MNColor()
        mix = g.Mix(
            factor_float=OutlineMask(threshold=threshold, thickness=thickness),
            a_color=group.o.color,
            b_color=(0.0, 0.0, 0.0, 1.0),
            data_type="RGBA",
            clamp_factor=True,
        )
        emission_1 = s.Emission(
            color=s.MenuSwitch.color(
                outline, {"Outline": mix.o.result_color, "None": group.o.color}
            ).o.output
        )

        emission_1 >> emission

        outline.default_value = "Outline"


ASSET = Flat
