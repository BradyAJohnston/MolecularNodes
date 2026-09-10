# Node group '.MN_utils_aa_atom_pos' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    CustomGeometryGroup,
    FloatSocket,
    IntegerSocket,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputInteger
from .utils_group_field_at_selection import Utils_group_field_at_selection


class MN_utils_aa_atom_pos(CustomGeometryGroup):
    """
    .MN_utils_aa_atom_pos

    Parameters
    ----------
    atom_name : InputInteger
        atom_name

    Inputs
    ------
    i.atom_name : IntegerSocket
        atom_name

    Outputs
    -------
    o.position : VectorSocket
        Position
    o.group_index : IntegerSocket
        Group Index
    o.b_factor : FloatSocket
        b_factor
    o.integer : IntegerSocket
        Integer
    """

    _name = ".MN_utils_aa_atom_pos"
    _tree_properties = {"node_tool_idname": "geometry._mn_utils_aa_atom_pos"}

    class _Inputs(SocketAccessor):
        atom_name: IntegerSocket

    class _Outputs(SocketAccessor):
        position: VectorSocket
        """Position"""
        group_index: IntegerSocket
        """Group Index"""
        b_factor: FloatSocket
        integer: IntegerSocket
        """Integer"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        atom_name: InputInteger = 5,
    ):
        super().__init__(**{"atom_name": atom_name})

    def _build_group(self, tree):
        atom_name = tree.inputs.integer("atom_name", 5)
        position = tree.outputs.vector("Position")
        group_index = tree.outputs.integer("Group Index")
        b_factor = tree.outputs.float("b_factor")
        integer = tree.outputs.integer("Integer")

        named_attribute = g.NamedAttribute.integer("atom_name")
        accumulate_field = g.AccumulateField.point.integer(
            g.Compare.integer.equal(named_attribute.o.attribute, 1)
        )
        group = Utils_group_field_at_selection(
            selection=g.Compare.integer.equal(named_attribute.o.attribute, atom_name),
            group_index=accumulate_field.o.leading,
            float=g.NamedAttribute.float("b_factor").o.attribute,
            vector=g.Position(),
            integer=g.EdgesOfVertex(vertex_index=g.Index()).o.total,
        )
        with g.Frame("If atom_name is 0, return midpoint of backbone N and C"):
            position_1 = g.Position()
            group_1 = Utils_group_field_at_selection(
                selection=g.Compare.integer.equal(named_attribute.o.attribute, 1),
                group_index=accumulate_field.o.leading,
                vector=position_1,
            )
            group_2 = Utils_group_field_at_selection(
                selection=g.Compare.integer.equal(named_attribute.o.attribute, 3),
                group_index=accumulate_field.o.leading,
                vector=position_1,
            )
            mix = g.Mix(
                a_vector=group_1.o.vector,
                b_vector=group_2.o.vector,
                factor_float=0.5,
                data_type="VECTOR",
                clamp_factor=True,
            )
            (
                g.Compare.integer.not_equal(atom_name, 0).o.result.switch.vector(
                    mix.o.result_vector, group.o.vector
                )
                >> position
            )

        accumulate_field >> group_index
        group.o.float >> b_factor
        group.o.integer >> integer
