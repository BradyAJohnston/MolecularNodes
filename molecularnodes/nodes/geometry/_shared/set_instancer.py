# Node group 'Set Instancer' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    BooleanSocket,
    CustomGeometryGroup,
    GeometrySocket,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputGeometry


class SetInstancer(CustomGeometryGroup):
    """
    Set Instancer

    Parameters
    ----------
    geometry : InputGeometry
        Atomic geometry that contains vertices and edges
    is_instanced : InputBoolean
        is_instanced

    Inputs
    ------
    i.geometry : GeometrySocket
        Atomic geometry that contains vertices and edges
    i.is_instanced : BooleanSocket
        is_instanced

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "Set Instancer"
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Atomic geometry that contains vertices and edges"""
        is_instanced: BooleanSocket

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
        is_instanced: InputBoolean = True,
    ):
        super().__init__(**{"Geometry": geometry, "is_instanced": is_instanced})

    def _build_group(self, tree):
        geometry = tree.inputs.geometry(
            "Geometry", description="Atomic geometry that contains vertices and edges"
        )
        is_instanced = tree.inputs.boolean("is_instanced", True)
        geometry_1 = tree.outputs.geometry("Geometry")

        (
            geometry
            >> g.StoreNamedAttribute.point.boolean(
                name="is_instanced", value=is_instanced
            )
            >> geometry_1
        )
