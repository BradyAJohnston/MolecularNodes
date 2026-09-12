# Node group ".Sample Nucleic Base Values" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    BooleanSocket,
    ColorSocket,
    CustomGeometryGroup,
    IntegerSocket,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputInteger
from ..color import Color
from ..residue_mask import ResidueMask
from ..select_nucleic_type import SelectNucleicType


class SampleNucleicBaseValues(CustomGeometryGroup):
    """
    .Sample Nucleic Base Values

    Parameters
    ----------
    input : InputInteger
        Input

    Inputs
    ------
    i.input : IntegerSocket
        Input

    Outputs
    -------
    o.base_valid : BooleanSocket
        base_valid
    o.base_pivot : VectorSocket
        base_pivot'
    o.base_z : VectorSocket
        base_Z
    o.base_y : VectorSocket
        base_Y
    o.base_position : VectorSocket
        base_position
    o.base_color : ColorSocket
        Base Color
    """

    _name = ".Sample Nucleic Base Values"
    _tree_properties = {"node_tool_idname": "geometry._sample_nucleic_base_values"}

    class _Inputs(SocketAccessor):
        input: IntegerSocket
        """Input"""

    class _Outputs(SocketAccessor):
        base_valid: BooleanSocket
        base_pivot: VectorSocket
        """base_pivot'"""
        base_z: VectorSocket
        """base_Z"""
        base_y: VectorSocket
        """base_Y"""
        base_position: VectorSocket
        base_color: ColorSocket
        """Base Color"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        input: InputInteger = 0,
    ):
        super().__init__(**{"Input": input})

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        _input = tree.inputs.integer("Input", 0)
        base_valid = tree.outputs.boolean("base_valid")
        base_pivot = tree.outputs.vector("base_pivot'")
        base_z = tree.outputs.vector("base_Z")
        base_y = tree.outputs.vector("base_Y")
        base_position = tree.outputs.vector("base_position")
        base_color = tree.outputs.color("Base Color", (0.8, 0.8, 0.8, 1.0))

        with g.Frame("Sample relevant base positions for orientations"):
            group = ResidueMask(atom_name=61)
            Color(index=ResidueMask(atom_name=67).o.index) >> base_color
            mix = g.Mix(
                a_vector=ResidueMask(atom_name=55).o.position,
                b_vector=ResidueMask(atom_name=57).o.position,
                data_type="VECTOR",
                clamp_factor=True,
            )
            group_1 = SelectNucleicType()
            group_2 = ResidueMask(
                atom_name=group_1.o.is_pyrimidine.switch.integer(65, 68)
            )
            group_3 = ResidueMask(
                atom_name=group_1.o.is_pyrimidine.switch.integer(62, 64)
            )
            (group.o.is_valid & group_2.o.is_valid & group_3.o.is_valid) >> base_valid
            group_3.o.position - group.o.position >> base_z
            group_2.o.position - group_3.o.position >> base_y

        group.o.position >> base_pivot
        mix.o.result_vector >> base_position
