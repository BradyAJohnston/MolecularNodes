# Node-group asset "Backbone O" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy.builder import (
    AssetGeometryGroup,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputMenu
from ._shared.backbone_position import BackbonePosition


class BackboneO(AssetGeometryGroup):
    """
    Backbone O

    Parameters
    ----------
    method : InputMenu | Literal["Read", "Compute"]
        Method

    Inputs
    ------
    i.method : MenuSocket
        Method

    Outputs
    -------
    o.o : VectorSocket
        O
    """

    _name = "Backbone O"
    _asset_name = "Backbone O"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        method: MenuSocket
        """Method"""

    class _Outputs(SocketAccessor):
        o: VectorSocket
        """O"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        method: InputMenu | Literal["Read", "Compute"] = "Compute",
    ):
        super().__init__(**{"Method": method})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        method = tree.inputs.menu("Method", expanded=True, optional_label=True)
        o = tree.outputs.vector("O")

        BackbonePosition(method=method, menu="backbone_O") >> o

        method.default_value = "Compute"


ASSET = BackboneO

ASSET_METADATA = {
    "catalog_id": "c5d9cc4a-2d12-48aa-838f-8114381e9e69",
}
