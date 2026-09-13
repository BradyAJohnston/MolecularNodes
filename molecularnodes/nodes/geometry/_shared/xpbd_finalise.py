# Node group "XPBD Finalise" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
)
from nodebpy.types import InputBoolean, InputFloat, InputGeometry


class XPBDFinalise(CustomGeometryGroup):
    """
    XPBD Finalise

    Parameters
    ----------
    geometry : InputGeometry
        Geometry
    selection : InputBoolean
        Selection
    deltat : InputFloat
        deltaT

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry
    i.selection : BooleanSocket
        Selection
    i.deltat : FloatSocket
        deltaT

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "XPBD Finalise"
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""
        selection: BooleanSocket
        """Selection"""
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
        deltat: InputFloat = 0.0,
    ):
        super().__init__(
            **{"Geometry": geometry, "Selection": selection, "deltaT": deltat}
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry("Geometry")
        selection = tree.inputs.boolean("Selection", True, hide_value=True)
        deltat = tree.inputs.float(
            "deltaT", 0.0, structure_type="SINGLE", force_non_field=True
        )
        geometry_1 = tree.outputs.geometry("Geometry")

        (
            geometry
            >> g.StoreNamedAttribute.point.vector(
                selection=selection,
                name="velocity",
                value=(
                    g.Position().o.position - g.NamedAttribute.vector("p_i").o.attribute
                )
                / deltat,
            )
            >> geometry_1
        )
