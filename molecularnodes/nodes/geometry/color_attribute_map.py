# Node-group asset "Color Attribute Map" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    ColorSocket,
    FloatSocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
    StringSocket,
)
from nodebpy.types import InputBoolean, InputColor, InputFloat, InputMenu, InputString
from .color_mix_intermediate import ColorMixIntermediate


class ColorAttributeMap(AssetGeometryGroup):
    """
    Color Attribute Map

    Parameters
    ----------
    color_space : InputMenu | Literal["Linear", "OKLab"]
        Color Space
    name : InputString
        Name of the attribute to map colors to
    min : InputFloat
        Value for the attribute to be the minimum color
    max : InputFloat
        Value for the attribute to be the maxium color
    socket_5 : InputBoolean
        Wheter to interpolate through the 'Mid' color.
    a : InputColor
        Color mapped to the minimum value of the attribute
    socket_7 : InputColor
        Color mapped to the middle value of the attribute
    b : InputColor
        Color mapped to the maximum value of the attribute

    Inputs
    ------
    i.color_space : MenuSocket
        Color Space
    i.name : StringSocket
        Name of the attribute to map colors to
    i.min : FloatSocket
        Value for the attribute to be the minimum color
    i.max : FloatSocket
        Value for the attribute to be the maxium color
    i.socket_5 : BooleanSocket
        Wheter to interpolate through the 'Mid' color.
    i.a : ColorSocket
        Color mapped to the minimum value of the attribute
    i.socket_7 : ColorSocket
        Color mapped to the middle value of the attribute
    i.b : ColorSocket
        Color mapped to the maximum value of the attribute

    Outputs
    -------
    o.color : ColorSocket
        The mapped color value based on the attribute.
    """

    _name = "Color Attribute Map"
    _asset_name = "Color Attribute Map"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "COLOR"
    _tree_properties = {"node_tool_idname": "geometry.color_attribute_map"}

    class _Inputs(SocketAccessor):
        color_space: MenuSocket
        """Color Space"""
        name: StringSocket
        """Name of the attribute to map colors to"""
        min: FloatSocket
        """Value for the attribute to be the minimum color"""
        max: FloatSocket
        """Value for the attribute to be the maxium color"""
        socket_5: BooleanSocket
        """Wheter to interpolate through the 'Mid' color."""
        a: ColorSocket
        """Color mapped to the minimum value of the attribute"""
        socket_7: ColorSocket
        """Color mapped to the middle value of the attribute"""
        b: ColorSocket
        """Color mapped to the maximum value of the attribute"""

    class _Outputs(SocketAccessor):
        color: ColorSocket
        """The mapped color value based on the attribute."""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        color_space: InputMenu | Literal["Linear", "OKLab"] = "Linear",
        name: InputString = "b_factor",
        min: InputFloat = 0.0,
        max: InputFloat = 150.0,
        socket_5: InputBoolean = True,
        a: InputColor = None,
        socket_7: InputColor = None,
        b: InputColor = None,
    ):
        super().__init__(
            **{
                "Color Space": color_space,
                "Name": name,
                "Min": min,
                "Max": max,
                "A": a,
                "B": b,
            },
            _named_links=[("Intermediate", socket_5), ("Intermediate", socket_7)],
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        color_space = tree.inputs.menu("Color Space", optional_label=True)
        name = tree.inputs.string(
            "Name",
            "b_factor",
            description="Name of the attribute to map colors to",
            optional_label=True,
        )
        min = tree.inputs.float(
            "Min",
            0.0,
            description="Value for the attribute to be the minimum color",
            min_value=-10_000.0,
            max_value=10_000.0,
        )
        max = tree.inputs.float(
            "Max",
            150.0,
            description="Value for the attribute to be the maxium color",
            min_value=-10_000.0,
            max_value=10_000.0,
        )
        with tree.inputs.panel("Color"):
            intermediate = tree.inputs.boolean(
                "Intermediate",
                True,
                description="Wheter to interpolate through the 'Mid' color.",
            )
            a = tree.inputs.color(
                "A",
                (0.0769495, 0.4785124, 0.5, 1.0),
                description="Color mapped to the minimum value of the attribute",
            )
            intermediate_1 = tree.inputs.color(
                "Intermediate",
                (0.5, 0.5, 0.5, 1.0),
                description="Color mapped to the middle value of the attribute",
            )
            b = tree.inputs.color(
                "B",
                (0.5, 0.1594808, 0.05802507, 1.0),
                description="Color mapped to the maximum value of the attribute",
            )
        color = tree.outputs.color(
            "Color",
            (0.0, 0.0, 0.0, 0.0),
            description="The mapped color value based on the attribute.",
        )

        (
            ColorMixIntermediate(
                factor=g.NamedAttribute.float(name).o.attribute.map_range(min, max),
                menu=color_space,
                socket_2=intermediate,
                a=a,
                socket_4=intermediate_1,
                b=b,
            )
            >> color
        )

        color_space.default_value = "Linear"


ASSET = ColorAttributeMap

ASSET_METADATA = {
    "catalog_id": "d3f975df-8408-4972-a669-8187a57e01d0",
}
