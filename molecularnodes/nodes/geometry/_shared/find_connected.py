# Node group 'Find Connected' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from nodebpy import geometry as g
from nodebpy.builder import (
    BooleanSocket,
    CustomGeometryGroup,
    IntegerSocket,
    MenuSocket,
    SocketAccessor,
)
from nodebpy.types import InputInteger, InputMenu


class FindConnected(CustomGeometryGroup):
    """
    Find Connected

    Parameters
    ----------
    value : InputInteger
        Value
    match : InputInteger
        Match
    distance : InputInteger
        Distance
    method : InputMenu | Literal["Any", "Exact"]
        Method

    Inputs
    ------
    i.value : IntegerSocket
        Value
    i.match : IntegerSocket
        Match
    i.distance : IntegerSocket
        Distance
    i.method : MenuSocket
        Method

    Outputs
    -------
    o.is_valid : BooleanSocket
        Is Valid
    o.length : IntegerSocket
        Length
    o.index : IntegerSocket
        Index
    """

    _name = "Find Connected"

    class _Inputs(SocketAccessor):
        value: IntegerSocket
        """Value"""
        match: IntegerSocket
        """Match"""
        distance: IntegerSocket
        """Distance"""
        method: MenuSocket
        """Method"""

    class _Outputs(SocketAccessor):
        is_valid: BooleanSocket
        """Is Valid"""
        length: IntegerSocket
        """Length"""
        index: IntegerSocket
        """Index"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        value: InputInteger = 0,
        match: InputInteger = 2,
        distance: InputInteger = 2,
        method: InputMenu | Literal["Any", "Exact"] = "Any",
    ):
        super().__init__(
            **{"Value": value, "Match": match, "Distance": distance, "Method": method}
        )

    def _build_group(self, tree):
        value = tree.inputs.integer("Value", 0, hide_value=True)
        match = tree.inputs.integer("Match", 2)
        distance = tree.inputs.integer("Distance", 2, min_value=1, max_value=3)
        method = tree.inputs.menu("Method", optional_label=True)
        is_valid = tree.outputs.boolean("Is Valid")
        length = tree.outputs.integer("Length")
        index = tree.outputs.integer("Index")

        shortest_edge_paths = g.ShortestEdgePaths(
            end_vertex=g.Compare.integer.equal(value, match)
        )
        compare = g.Compare.integer.equal(shortest_edge_paths.o.total_cost, 2)
        compare_1 = g.Compare.integer.equal(shortest_edge_paths.o.total_cost, 3)
        evaluate_at_index = shortest_edge_paths.o.next_vertex_index.point.at(
            shortest_edge_paths.o.next_vertex_index
        )
        evaluate_at_index_1 = shortest_edge_paths.o.next_vertex_index.point.at(
            evaluate_at_index
        )
        float_to_integer = g.FloatToInteger(float=shortest_edge_paths.o.total_cost)
        switch = g.Compare.integer.equal(
            shortest_edge_paths.o.total_cost, 1
        ).o.result.switch.integer(-1, shortest_edge_paths.o.next_vertex_index)
        switch_1 = compare.o.result.switch.integer(switch, evaluate_at_index)
        index_switch = g.IndexSwitch.integer(
            distance,
            (
                0,
                switch,
                switch_1,
                compare_1.o.result.switch.integer(switch_1, evaluate_at_index_1),
            ),
        )
        index_switch_1 = g.IndexSwitch.integer(
            distance,
            (
                0,
                switch,
                compare.o.result.switch.integer(-1, evaluate_at_index),
                compare_1.o.result.switch.integer(-1, evaluate_at_index_1),
            ),
        )
        menu_switch = g.MenuSwitch.integer(
            method, {"Any": index_switch, "Exact": index_switch_1}
        )
        g.Compare.integer.not_equal(menu_switch.o.output, -1) >> is_valid

        float_to_integer >> length
        menu_switch >> index

        method.default_value = "Any"
