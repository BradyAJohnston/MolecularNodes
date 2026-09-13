# Node-group asset "World to Angstrom" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputFloat
from ._shared.mn_world_scale import MN_world_scale


class WorldToAngstrom(AssetGeometryGroup):
    """
    World to Angstrom

    Parameters
    ----------
    world : InputFloat
        World

    Inputs
    ------
    i.world : FloatSocket
        World

    Outputs
    -------
    o.angstrom : FloatSocket
        Angstrom
    """

    _name = "World to Angstrom"
    _asset_name = "World to Angstrom"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {"node_tool_idname": "geometry.world_to_angstrom"}

    class _Inputs(SocketAccessor):
        world: FloatSocket
        """World"""

    class _Outputs(SocketAccessor):
        angstrom: FloatSocket
        """Angstrom"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        world: InputFloat = 0.5,
    ):
        super().__init__(**{"World": world})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        world = tree.inputs.float("World", 0.5, min_value=-10_000.0, max_value=10_000.0)
        angstrom = tree.outputs.float("Angstrom")

        world / MN_world_scale() >> angstrom


ASSET = WorldToAngstrom

ASSET_METADATA = {
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
