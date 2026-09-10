# Node-group asset 'Transform Mix' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    MatrixSocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputFloat, InputMatrix, InputMenu


class TransformMix(AssetGeometryGroup):
    """
    Mix between two transformations

    Parameters
    ----------
    a : InputMatrix
        Transform A to mix from at 0.0
    b : InputMatrix
        Transform B which will be mixed to at 1.0
    menu : InputMenu | Literal["Single", "Split"]
        Menu
    translation : InputFloat
        Amount to mix the Translation between A and B
    rotation : InputFloat
        Amount to mix the Rotation between A and B
    scale : InputFloat
        Amount to mix the Scale between A and B
    factor : InputFloat
        Factor

    Inputs
    ------
    i.a : MatrixSocket
        Transform A to mix from at 0.0
    i.b : MatrixSocket
        Transform B which will be mixed to at 1.0
    i.menu : MenuSocket
        Menu
    i.translation : FloatSocket
        Amount to mix the Translation between A and B
    i.rotation : FloatSocket
        Amount to mix the Rotation between A and B
    i.scale : FloatSocket
        Amount to mix the Scale between A and B
    i.factor : FloatSocket
        Factor

    Outputs
    -------
    o.transform : MatrixSocket
        The final mixed Transform
    """

    _name = "Transform Mix"
    _asset_name = "Transform Mix"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {
        "description": "Mix between two transformations",
        "node_tool_idname": "geometry.transform_mix",
    }

    class _Inputs(SocketAccessor):
        a: MatrixSocket
        """Transform A to mix from at 0.0"""
        b: MatrixSocket
        """Transform B which will be mixed to at 1.0"""
        menu: MenuSocket
        """Menu"""
        translation: FloatSocket
        """Amount to mix the Translation between A and B"""
        rotation: FloatSocket
        """Amount to mix the Rotation between A and B"""
        scale: FloatSocket
        """Amount to mix the Scale between A and B"""
        factor: FloatSocket
        """Factor"""

    class _Outputs(SocketAccessor):
        transform: MatrixSocket
        """The final mixed Transform"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        a: InputMatrix = None,
        b: InputMatrix = None,
        menu: InputMenu | Literal["Single", "Split"] = "Single",
        translation: InputFloat = 0.5,
        rotation: InputFloat = 0.5,
        scale: InputFloat = 0.5,
        factor: InputFloat = 0.5,
    ):
        super().__init__(
            **{
                "A": a,
                "B": b,
                "Menu": menu,
                "Translation": translation,
                "Rotation": rotation,
                "Scale": scale,
                "Factor": factor,
            }
        )

    def _build_group(self, tree):
        a = tree.inputs.matrix("A", description="Transform A to mix from at 0.0")
        b = tree.inputs.matrix(
            "B", description="Transform B which will be mixed to at 1.0"
        )
        menu = tree.inputs.menu("Menu", expanded=True, optional_label=True)
        translation = tree.inputs.float(
            "Translation",
            0.5,
            description="Amount to mix the Translation between A and B",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        rotation = tree.inputs.float(
            "Rotation",
            0.5,
            description="Amount to mix the Rotation between A and B",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        scale = tree.inputs.float(
            "Scale",
            0.5,
            description="Amount to mix the Scale between A and B",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        factor = tree.inputs.float(
            "Factor", 0.5, min_value=0.0, max_value=1.0, subtype="FACTOR"
        )
        transform = tree.outputs.matrix(
            "Transform", description="The final mixed Transform"
        )

        menu_switch = g.MenuSwitch.integer(menu, {"Single": 0, "Split": 1})
        mix = g.Mix(
            factor_float=g.IndexSwitch.float(menu_switch.o.output, (factor, rotation)),
            a_rotation=a.rotation,
            b_rotation=b.rotation,
            data_type="ROTATION",
            clamp_factor=True,
        )
        mix_1 = g.Mix(
            factor_float=g.IndexSwitch.float(
                menu_switch.o.output, (factor, translation)
            ),
            a_vector=a.translation,
            b_vector=b.translation,
            data_type="VECTOR",
            clamp_factor=True,
        )
        mix_2 = g.Mix(
            factor_float=g.IndexSwitch.float(menu_switch.o.output, (factor, scale)),
            a_vector=a.scale,
            b_vector=b.scale,
            data_type="VECTOR",
            clamp_factor=True,
        )
        combine_transform = g.CombineTransform(
            translation=mix_1.o.result_vector,
            rotation=mix.o.result_rotation,
            scale=mix_2.o.result_vector,
        )

        combine_transform >> transform

        menu.default_value = "Single"


ASSET = TransformMix

ASSET_METADATA = {
    "description": "Mix between two transformations",
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
