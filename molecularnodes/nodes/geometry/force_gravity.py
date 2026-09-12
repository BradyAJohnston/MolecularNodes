# Node-group asset "Force Gravity" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    AssetGeometryGroup,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputVector


class ForceGravity(AssetGeometryGroup):
    """
    Force Gravity

    Parameters
    ----------
    add : InputVector
        Add
    gravity : InputVector
        Gravity

    Inputs
    ------
    i.add : VectorSocket
        Add
    i.gravity : VectorSocket
        Gravity

    Outputs
    -------
    o.force : VectorSocket
        Force
    """

    _name = "Force Gravity"
    _asset_name = "Force Gravity"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        add: VectorSocket
        """Add"""
        gravity: VectorSocket
        """Gravity"""

    class _Outputs(SocketAccessor):
        force: VectorSocket
        """Force"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        add: InputVector = None,
        gravity: InputVector = None,
    ):
        super().__init__(**{"Add": add, "Gravity": gravity})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        add = tree.inputs.vector(
            "Add",
            (0.0, 0.0, 0.0),
            min_value=-10_000.0,
            max_value=10_000.0,
            hide_value=True,
        )
        gravity = tree.inputs.vector(
            "Gravity", (0.0, 0.0, -9.8), min_value=-10_000.0, max_value=10_000.0
        )
        force = tree.outputs.vector("Force")

        add + gravity >> force


ASSET = ForceGravity

ASSET_METADATA = {
    "catalog_id": "c2c958af-5095-4fc2-884d-709bba965fc4",
}
