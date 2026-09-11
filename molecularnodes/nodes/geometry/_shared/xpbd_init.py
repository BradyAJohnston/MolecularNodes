# Node group "XPBD Init" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    BooleanSocket,
    CustomGeometryGroup,
    FloatSocket,
    GeometrySocket,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputBoolean, InputFloat, InputGeometry, InputVector
from ..velocity import Velocity
from .inverse_mass import InverseMass


class XPBDInit(CustomGeometryGroup):
    """
    XPBD Init

    Parameters
    ----------
    geometry : InputGeometry
        Geometry
    selection : InputBoolean
        Selection
    force : InputVector
        Force
    drag : InputFloat
        Drag
    deltat : InputFloat
        deltaT

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry
    i.selection : BooleanSocket
        Selection
    i.force : VectorSocket
        Force
    i.drag : FloatSocket
        Drag
    i.deltat : FloatSocket
        deltaT

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "XPBD Init"
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""
        selection: BooleanSocket
        """Selection"""
        force: VectorSocket
        """Force"""
        drag: FloatSocket
        """Drag"""
        deltat: FloatSocket
        """deltaT"""

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
        selection: InputBoolean = True,
        force: InputVector = None,
        drag: InputFloat = 10.0,
        deltat: InputFloat = 1.0,
    ):
        super().__init__(
            **{
                "Geometry": geometry,
                "Selection": selection,
                "Force": force,
                "Drag": drag,
                "deltaT": deltat,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry("Geometry")
        selection = tree.inputs.boolean("Selection", True, hide_value=True)
        force = tree.inputs.vector(
            "Force", (0.0, 0.0, -9.8), min_value=-10_000.0, max_value=10_000.0
        )
        drag = tree.inputs.float("Drag", 10.0, min_value=-10_000.0, max_value=10_000.0)
        deltat = tree.inputs.float(
            "deltaT",
            1.0,
            min_value=-10_000.0,
            max_value=10_000.0,
            structure_type="SINGLE",
            force_non_field=True,
        )
        geometry_1 = tree.outputs.geometry("Geometry")

        with g.Frame("Drag Force"):
            vector_math = Velocity().o.velocity * (1.0 - drag * deltat).clamp()
        with g.Frame("New Forces"):
            vector_math_1 = force * deltat * InverseMass().o.w
        (
            geometry
            >> g.StoreNamedAttribute.point.vector(name="p_i", value=g.Position())
            >> g.StoreNamedAttribute.point.vector(
                name="velocity", value=vector_math + vector_math_1
            )
            >> g.SetPosition(selection=selection, offset=Velocity().o.velocity * deltat)
            >> geometry_1
        )
