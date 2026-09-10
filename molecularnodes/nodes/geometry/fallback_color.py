# Node-group asset 'Fallback Color' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    ColorSocket,
    PackageLibrary,
    SocketAccessor,
    StringSocket,
)
from nodebpy.types import InputColor, InputString


class FallbackColor(AssetGeometryGroup):
    """
    Fallback Color

    Parameters
    ----------
    name : InputString
        Name of the attribute to attempt to read from the geometry
    fallback : InputColor
        Value to use instead if the named attribute doesn't exist on the geometry

    Inputs
    ------
    i.name : StringSocket
        Name of the attribute to attempt to read from the geometry
    i.fallback : ColorSocket
        Value to use instead if the named attribute doesn't exist on the geometry

    Outputs
    -------
    o.color : ColorSocket
        The named attribute read from the geometry if it exists, or the fallback value if it doesn't
    """

    _name = "Fallback Color"
    _asset_name = "Fallback Color"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.fallback_color"}

    class _Inputs(SocketAccessor):
        name: StringSocket
        """Name of the attribute to attempt to read from the geometry"""
        fallback: ColorSocket
        """Value to use instead if the named attribute doesn't exist on the geometry"""

    class _Outputs(SocketAccessor):
        color: ColorSocket
        """The named attribute read from the geometry if it exists, or the fallback value if it doesn't"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        name: InputString = "Color",
        fallback: InputColor = None,
    ):
        super().__init__(**{"Name": name, "Fallback": fallback})

    def _build_group(self, tree):
        name = tree.inputs.string(
            "Name",
            "Color",
            description="Name of the attribute to attempt to read from the geometry",
            optional_label=True,
        )
        fallback = tree.inputs.color(
            "Fallback",
            (0.0716495, 0.2954248, 0.23565927, 1.0),
            description="Value to use instead if the named attribute doesn't exist on the geometry",
        )
        color = tree.outputs.color(
            "Color",
            (0.0, 0.0, 0.0, 1.0),
            description="The named attribute read from the geometry if it exists, or the fallback value if it doesn't",
        )

        named_attribute = g.NamedAttribute.color(name)
        (
            named_attribute.o.exists.switch.color(fallback, named_attribute.o.attribute)
            >> color
        )


ASSET = FallbackColor

ASSET_METADATA = {
    "catalog_id": "7ccb8802-a69f-483e-bf6e-4a47aaa9e940",
}
