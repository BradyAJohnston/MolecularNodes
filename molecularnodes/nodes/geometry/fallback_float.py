# Node-group asset 'Fallback Float' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    PackageLibrary,
    SocketAccessor,
    StringSocket,
)
from nodebpy.types import InputFloat, InputString


class FallbackFloat(AssetGeometryGroup):
    """
    Fallback Float

    Parameters
    ----------
    name : InputString
        Name of the attribute to attempt to read from the geometry
    fallback : InputFloat
        Value to use instead if the named attribute doesn't exist on the geometry

    Inputs
    ------
    i.name : StringSocket
        Name of the attribute to attempt to read from the geometry
    i.fallback : FloatSocket
        Value to use instead if the named attribute doesn't exist on the geometry

    Outputs
    -------
    o.value : FloatSocket
        The named attribute read from the geometry if it exists, or the fallback value if it doesn't
    """

    _name = "Fallback Float"
    _asset_name = "Fallback Float"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.fallback_float"}

    class _Inputs(SocketAccessor):
        name: StringSocket
        """Name of the attribute to attempt to read from the geometry"""
        fallback: FloatSocket
        """Value to use instead if the named attribute doesn't exist on the geometry"""

    class _Outputs(SocketAccessor):
        value: FloatSocket
        """The named attribute read from the geometry if it exists, or the fallback value if it doesn't"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        name: InputString = "",
        fallback: InputFloat = 0.0,
    ):
        super().__init__(**{"Name": name, "Fallback": fallback})

    def _build_group(self, tree):
        name = tree.inputs.string(
            "Name",
            "",
            description="Name of the attribute to attempt to read from the geometry",
            optional_label=True,
        )
        fallback = tree.inputs.float(
            "Fallback",
            0.0,
            description="Value to use instead if the named attribute doesn't exist on the geometry",
        )
        value = tree.outputs.float(
            "Value",
            description="The named attribute read from the geometry if it exists, or the fallback value if it doesn't",
        )

        named_attribute = g.NamedAttribute.float(name)
        (
            named_attribute.o.exists.switch.float(fallback, named_attribute.o.attribute)
            >> value
        )


ASSET = FallbackFloat

ASSET_METADATA = {
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
