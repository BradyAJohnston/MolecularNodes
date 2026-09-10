# Node-group asset 'Select Bonded' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputInteger
from .boolean_andor import BooleanAndOr


class SelectBonded(AssetGeometryGroup):
    """
    Select Bonded

    Parameters
    ----------
    selection : InputBoolean
        Selection of atoms to apply this node to
    depth : InputInteger
        Number of bonds to expand the selection by

    Inputs
    ------
    i.selection : BooleanSocket
        Selection of atoms to apply this node to
    i.depth : IntegerSocket
        Number of bonds to expand the selection by

    Outputs
    -------
    o.selection : BooleanSocket
        The calculated selection
    o.bonded : BooleanSocket
        Expanded Selection that excludes the original selection
    """

    _name = "Select Bonded"
    _asset_name = "Select Bonded"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"
    _tree_properties = {"node_tool_idname": "geometry.select_bonded"}

    class _Inputs(SocketAccessor):
        selection: BooleanSocket
        """Selection of atoms to apply this node to"""
        depth: IntegerSocket
        """Number of bonds to expand the selection by"""

    class _Outputs(SocketAccessor):
        selection: BooleanSocket
        """The calculated selection"""
        bonded: BooleanSocket
        """Expanded Selection that excludes the original selection"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        selection: InputBoolean = False,
        depth: InputInteger = 1,
    ):
        super().__init__(**{"Selection": selection, "Depth": depth})

    def _build_group(self, tree):
        selection = tree.inputs.boolean(
            "Selection",
            False,
            description="Selection of atoms to apply this node to",
            hide_value=True,
        )
        depth = tree.inputs.integer(
            "Depth",
            1,
            description="Number of bonds to expand the selection by",
            min_value=0,
        )
        selection_1 = tree.outputs.boolean(
            "Selection", description="The calculated selection"
        )
        bonded = tree.outputs.boolean(
            "Bonded",
            description="Expanded Selection that excludes the original selection",
        )

        shortest_edge_paths = g.ShortestEdgePaths(end_vertex=selection)
        compare = g.Compare.integer.less_equal(shortest_edge_paths.o.total_cost, depth)
        compare_1 = g.Compare.integer.greater_than(shortest_edge_paths.o.total_cost, 0)
        (compare.o.result & compare_1) >> bonded
        BooleanAndOr(and_=compare, or_=selection, boolean=compare_1) >> selection_1


ASSET = SelectBonded

ASSET_METADATA = {
    "catalog_id": "bd1f205b-fea5-4700-b2c2-754f3321e969",
}
