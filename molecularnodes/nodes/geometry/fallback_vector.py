# Node-group asset "Fallback Vector" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    PackageLibrary,
    SocketAccessor,
    StringSocket,
    VectorSocket,
)
from nodebpy.types import InputString, InputVector


class FallbackVector(AssetGeometryGroup):
    """
    Fallback Vector

    Parameters
    ----------
    name : InputString
        Name of the attribute to attempt to read from the geometry
    fallback : InputVector
        Value to use instead if the named attribute doesn't exist on the geometry

    Inputs
    ------
    i.name : StringSocket
        Name of the attribute to attempt to read from the geometry
    i.fallback : VectorSocket
        Value to use instead if the named attribute doesn't exist on the geometry

    Outputs
    -------
    o.output : VectorSocket
        The named attribute read from the geometry if it exists, or the fallback value if it doesn't
    """

    _name = "Fallback Vector"
    _asset_name = "Fallback Vector"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.fallback_vector"}

    class _Inputs(SocketAccessor):
        name: StringSocket
        """Name of the attribute to attempt to read from the geometry"""
        fallback: VectorSocket
        """Value to use instead if the named attribute doesn't exist on the geometry"""

    class _Outputs(SocketAccessor):
        output: VectorSocket
        """The named attribute read from the geometry if it exists, or the fallback value if it doesn't"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        name: InputString = "",
        fallback: InputVector = None,
    ):
        super().__init__(**{"Name": name, "Fallback": fallback})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        name = tree.inputs.string(
            "Name",
            "",
            description="Name of the attribute to attempt to read from the geometry",
            optional_label=True,
        )
        fallback = tree.inputs.vector(
            "Fallback",
            (0.0, 0.0, 0.0),
            description="Value to use instead if the named attribute doesn't exist on the geometry",
        )
        output = tree.outputs.vector(
            "Output",
            description="The named attribute read from the geometry if it exists, or the fallback value if it doesn't",
        )

        named_attribute = g.NamedAttribute.vector(name)
        (
            named_attribute.o.exists.switch.vector(
                fallback, named_attribute.o.attribute
            )
            >> output
        )


ASSET = FallbackVector

ASSET_METADATA = {
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
