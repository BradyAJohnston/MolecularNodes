# Node group "XPBD Solve Hook" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from .lerp_position import LerpPosition


class XPBDSolveHook(CustomGeometryGroup):
    """
    XPBD Solve Hook

    Parameters
    ----------
    geometry : InputGeometry
        Geometry
    selection : InputBoolean
        Selection
    target : InputVector
        Target
    decay : InputFloat
        Decay
    deltat : InputFloat
        deltaT

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry
    i.selection : BooleanSocket
        Selection
    i.target : VectorSocket
        Target
    i.decay : FloatSocket
        Decay
    i.deltat : FloatSocket
        deltaT

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "XPBD Solve Hook"
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""
        selection: BooleanSocket
        """Selection"""
        target: VectorSocket
        """Target"""
        decay: FloatSocket
        """Decay"""
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
        selection: InputBoolean = False,
        target: InputVector = None,
        decay: InputFloat = 0.5,
        deltat: InputFloat = 0.5,
    ):
        super().__init__(
            **{
                "Geometry": geometry,
                "Selection": selection,
                "Target": target,
                "Decay": decay,
                "deltaT": deltat,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry("Geometry")
        selection = tree.inputs.boolean("Selection", False, hide_value=True)
        target = tree.inputs.vector(
            "Target", (0.0, 0.0, 0.0), min_value=-10_000.0, max_value=10_000.0
        )
        decay = tree.inputs.float("Decay", 0.5, min_value=-10_000.0, max_value=10_000.0)
        deltat = tree.inputs.float(
            "deltaT",
            0.5,
            min_value=-10_000.0,
            max_value=10_000.0,
            structure_type="SINGLE",
            force_non_field=True,
        )
        geometry_1 = tree.outputs.geometry("Geometry")

        group = LerpPosition(b=target, deltat=deltat, decay=decay)
        (
            geometry
            >> g.SetPosition(selection=selection, position=group.o.position)
            >> geometry_1
        )
        _switch = (group.o.length < 0.01).switch.vector(group.o.position, target)
