# Node group ".MN_bs_smooth" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    CustomGeometryGroup,
    FloatSocket,
    GeometrySocket,
    IntegerSocket,
    SocketAccessor,
)
from nodebpy.types import InputFloat, InputGeometry, InputInteger
from ..expand_boolean import ExpandBoolean
from .mn_select_sec_struct import MN_select_sec_struct


class MN_bs_smooth(CustomGeometryGroup):
    """
    .MN_bs_smooth

    Parameters
    ----------
    geometry : InputGeometry
        Geometry
    factor : InputFloat
        Factor
    iterations : InputInteger
        Iterations

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry
    i.factor : FloatSocket
        Factor
    i.iterations : IntegerSocket
        Iterations

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = ".MN_bs_smooth"
    _color_tag = "GEOMETRY"
    _tree_properties = {"node_tool_idname": "geometry._mn_bs_smooth"}

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""
        factor: FloatSocket
        """Factor"""
        iterations: IntegerSocket
        """Iterations"""

    class _Outputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        geometry: InputGeometry = None,
        factor: InputFloat = 1.0,
        iterations: InputInteger = 2,
    ):
        super().__init__(
            **{"Geometry": geometry, "Factor": factor, "Iterations": iterations}
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry("Geometry")
        factor = tree.inputs.float(
            "Factor", 1.0, min_value=0.0, max_value=1.0, subtype="FACTOR"
        )
        iterations = tree.inputs.integer("Iterations", 2, min_value=0)
        geometry_1 = tree.outputs.geometry("Geometry")

        group = MN_select_sec_struct()
        position = g.Position()
        blur_attribute = g.BlurAttribute.vector(
            position, iterations, ExpandBoolean(boolean=group.o.is_structured, expand=1)
        )
        mix = g.Mix(
            factor_float=factor,
            a_vector=position,
            b_vector=blur_attribute,
            data_type="VECTOR",
            clamp_factor=True,
        )
        (
            geometry
            >> g.SetPosition(selection=group.o.is_sheet, position=mix.o.result_vector)
            >> geometry_1
        )
