# Node-group asset "Fallback Boolean" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    PackageLibrary,
    SocketAccessor,
    StringSocket,
)
from nodebpy.types import InputBoolean, InputString


class FallbackBoolean(AssetGeometryGroup):
    """
    Computes the boolean field if the given attribute doesn't exist. If it doesn't exist it just uses the attribute instead

    Parameters
    ----------
    name : InputString
        Name of the attribute to attempt to read from the geometry
    fallback : InputBoolean
        Value to use instead if the named attribute doesn't exist on the geometry

    Inputs
    ------
    i.name : StringSocket
        Name of the attribute to attempt to read from the geometry
    i.fallback : BooleanSocket
        Value to use instead if the named attribute doesn't exist on the geometry

    Outputs
    -------
    o.boolean : BooleanSocket
        The named attribute read from the geometry if it exists, or the fallback value if it doesn't
    """

    _name = "Fallback Boolean"
    _asset_name = "Fallback Boolean"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {
        "description": "Computes the boolean field if the given attribute doesn't exist. If it doesn't exist it just uses the attribute instead",
        "node_tool_idname": "geometry.fallback_boolean",
    }

    class _Inputs(SocketAccessor):
        name: StringSocket
        """Name of the attribute to attempt to read from the geometry"""
        fallback: BooleanSocket
        """Value to use instead if the named attribute doesn't exist on the geometry"""

    class _Outputs(SocketAccessor):
        boolean: BooleanSocket
        """The named attribute read from the geometry if it exists, or the fallback value if it doesn't"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        name: InputString = "",
        fallback: InputBoolean = False,
    ):
        super().__init__(**{"Name": name, "Fallback": fallback})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        name = tree.inputs.string(
            "Name",
            "",
            description="Name of the attribute to attempt to read from the geometry",
            optional_label=True,
        )
        fallback = tree.inputs.boolean(
            "Fallback",
            False,
            description="Value to use instead if the named attribute doesn't exist on the geometry",
        )
        boolean = tree.outputs.boolean(
            "Boolean",
            description="The named attribute read from the geometry if it exists, or the fallback value if it doesn't",
        )

        named_attribute = g.NamedAttribute.boolean(name)
        (
            named_attribute.o.exists.switch.boolean(
                fallback, named_attribute.o.attribute
            )
            >> boolean
        )


ASSET = FallbackBoolean

ASSET_METADATA = {
    "description": "Computes the boolean field if the given attribute doesn't exist. If it doesn't exist it just uses the attribute instead",
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
