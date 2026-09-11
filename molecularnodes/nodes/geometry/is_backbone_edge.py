# Node-group asset "Is Backbone Edge" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
)
from .chain_id import ChainID
from .edge_group_id import EdgeGroupID
from .is_alpha_carbon import IsAlphaCarbon
from .residue_id import ResidueID


class IsBackboneEdge(AssetGeometryGroup):
    """
    Is Backbone Edge

    Outputs
    -------
    o.is_backbone_edge : BooleanSocket
        Is Backbone Edge
    """

    _name = "Is Backbone Edge"
    _asset_name = "Is Backbone Edge"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        pass

    class _Outputs(SocketAccessor):
        is_backbone_edge: BooleanSocket
        """Is Backbone Edge"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(self):
        super().__init__()

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        is_backbone_edge = tree.outputs.boolean(
            "Is Backbone Edge", attribute_domain="EDGE"
        )

        boolean_math = (
            g.Compare.integer.equal(
                EdgeGroupID(group_id=ResidueID()).o.difference, 1
            ).o.result
            & EdgeGroupID(group_id=ChainID()).o.is_equal
        )
        (boolean_math & IsAlphaCarbon().o.selection.edge.evaluate()) >> is_backbone_edge


ASSET = IsBackboneEdge

ASSET_METADATA = {
    "catalog_id": "bd1f205b-fea5-4700-b2c2-754f3321e969",
}
