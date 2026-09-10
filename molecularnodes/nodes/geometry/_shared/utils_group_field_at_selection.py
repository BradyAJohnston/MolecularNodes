# Node group '.utils_group_field_at_selection' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    BooleanSocket,
    ColorSocket,
    CustomGeometryGroup,
    FloatSocket,
    IntegerSocket,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import (
    InputBoolean,
    InputColor,
    InputFloat,
    InputInteger,
    InputVector,
)


class Utils_group_field_at_selection(CustomGeometryGroup):
    """
    .utils_group_field_at_selection

    Parameters
    ----------
    selection : InputBoolean
        Selection of atoms to apply this node to
    group_index : InputInteger
        Group Index
    float : InputFloat
        Float
    vector : InputVector
        Vector
    boolean : InputBoolean
        Boolean
    color : InputColor
        Color
    integer : InputInteger
        Integer

    Inputs
    ------
    i.selection : BooleanSocket
        Selection of atoms to apply this node to
    i.group_index : IntegerSocket
        Group Index
    i.float : FloatSocket
        Float
    i.vector : VectorSocket
        Vector
    i.boolean : BooleanSocket
        Boolean
    i.color : ColorSocket
        Color
    i.integer : IntegerSocket
        Integer

    Outputs
    -------
    o.group_index : IntegerSocket
        Group Index
    o.float : FloatSocket
        Float
    o.vector : VectorSocket
        Vector
    o.boolean : BooleanSocket
        Boolean
    o.color : ColorSocket
        Color
    o.integer : IntegerSocket
        Integer
    """

    _name = ".utils_group_field_at_selection"
    _tree_properties = {"node_tool_idname": "geometry._utils_group_field_at_selection"}

    class _Inputs(SocketAccessor):
        selection: BooleanSocket
        """Selection of atoms to apply this node to"""
        group_index: IntegerSocket
        """Group Index"""
        float: FloatSocket
        """Float"""
        vector: VectorSocket
        """Vector"""
        boolean: BooleanSocket
        """Boolean"""
        color: ColorSocket
        """Color"""
        integer: IntegerSocket
        """Integer"""

    class _Outputs(SocketAccessor):
        group_index: IntegerSocket
        """Group Index"""
        float: FloatSocket
        """Float"""
        vector: VectorSocket
        """Vector"""
        boolean: BooleanSocket
        """Boolean"""
        color: ColorSocket
        """Color"""
        integer: IntegerSocket
        """Integer"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        selection: InputBoolean = False,
        group_index: InputInteger = 0,
        float: InputFloat = 0.0,
        vector: InputVector = None,
        boolean: InputBoolean = False,
        color: InputColor = None,
        integer: InputInteger = 0,
    ):
        super().__init__(
            **{
                "Selection": selection,
                "Group Index": group_index,
                "Float": float,
                "Vector": vector,
                "Boolean": boolean,
                "Color": color,
                "Integer": integer,
            }
        )

    def _build_group(self, tree):
        selection = tree.inputs.boolean(
            "Selection", False, description="Selection of atoms to apply this node to"
        )
        group_index = tree.inputs.integer("Group Index", 0)
        float = tree.inputs.float("Float", 0.0, hide_value=True)
        vector = tree.inputs.vector("Vector", (0.0, 0.0, 0.0), hide_value=True)
        boolean = tree.inputs.boolean("Boolean", False, hide_value=True)
        color = tree.inputs.color("Color", (0.0, 0.0, 0.0, 0.0), hide_value=True)
        integer = tree.inputs.integer("Integer", 0, hide_value=True)
        group_index_1 = tree.outputs.integer("Group Index")
        float_1 = tree.outputs.float("Float")
        vector_1 = tree.outputs.vector("Vector")
        boolean_1 = tree.outputs.boolean("Boolean")
        color_1 = tree.outputs.color("Color", (0.0, 0.0, 0.0, 0.0))
        integer_1 = tree.outputs.integer("Integer")

        accumulate_field = selection.switch.integer(true=g.Index()).point.total(
            group_index
        )
        float.point.at(accumulate_field) >> float_1
        vector.point.at(accumulate_field) >> vector_1
        boolean.point.at(accumulate_field) >> boolean_1
        color.point.at(accumulate_field) >> color_1
        integer.point.at(accumulate_field) >> integer_1

        accumulate_field >> group_index_1
