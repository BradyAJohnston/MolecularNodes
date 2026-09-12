# Node-group asset "Angstrom to World" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from ._shared.unit_convert import UnitConvert


class AngstromToWorld(AssetGeometryGroup):
    """
    Angstrom to World

    Parameters
    ----------
    angstrom : InputFloat
        Angstrom

    Inputs
    ------
    i.angstrom : FloatSocket
        Angstrom

    Outputs
    -------
    o.world : FloatSocket
        World
    """

    _name = "Angstrom to World"
    _asset_name = "Angstrom to World"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "CONVERTER"
    _tree_properties = {"node_tool_idname": "geometry.angstrom_to_world"}

    class _Inputs(SocketAccessor):
        angstrom: FloatSocket
        """Angstrom"""

    class _Outputs(SocketAccessor):
        world: FloatSocket
        """World"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        angstrom: InputFloat = 3.0,
    ):
        super().__init__(**{"Angstrom": angstrom})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        angstrom = tree.inputs.float(
            "Angstrom", 3.0, min_value=-10_000.0, max_value=10_000.0
        )
        world = tree.outputs.float("World")

        UnitConvert(from_=angstrom) >> world


ASSET = AngstromToWorld

ASSET_METADATA = {
    "catalog_id": "b293127a-ef53-4981-b170-fce54963caa7",
}
