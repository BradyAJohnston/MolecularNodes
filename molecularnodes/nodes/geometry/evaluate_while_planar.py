# Node-group asset "Evaluate While Planar" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    ClosureSocket,
    GeometrySocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputClosure, InputGeometry
from .geometry_to_planar import GeometryToPlanar


class EvaluateWhilePlanar(AssetGeometryGroup):
    """
    Evaluate While Planar

    Parameters
    ----------
    geometry : InputGeometry
        Geometry to transform
    selection : InputBoolean
        The parts of the geometry that contribute to the planar calculation
    closure : InputClosure
        Closure

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry to transform
    i.selection : BooleanSocket
        The parts of the geometry that contribute to the planar calculation
    i.closure : ClosureSocket
        Closure

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "Evaluate While Planar"
    _asset_name = "Evaluate While Planar"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry to transform"""
        selection: BooleanSocket
        """The parts of the geometry that contribute to the planar calculation"""
        closure: ClosureSocket
        """Closure"""

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
        closure: InputClosure = None,
    ):
        super().__init__(
            **{"Geometry": geometry, "Selection": selection, "Closure": closure}
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry("Geometry", description="Geometry to transform")
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="The parts of the geometry that contribute to the planar calculation",
            hide_value=True,
        )
        closure = tree.inputs.closure("Closure")
        geometry_1 = tree.outputs.geometry("Geometry")

        geometry_to_planar = GeometryToPlanar(geometry=geometry, selection=selection)
        evaluate_closure = g.EvaluateClosure(closure)
        evaluate_closure.inputs.geometry(geometry_to_planar.o.geometry, "Geometry")
        geometry_2 = evaluate_closure.outputs.geometry("Geometry")
        (
            geometry_2.output
            >> g.TransformGeometry(
                transform=geometry_to_planar.o.transform.invert(), mode="Matrix"
            )
            >> geometry_1
        )


ASSET = EvaluateWhilePlanar

ASSET_METADATA = {
    "catalog_id": "a1e4128a-131f-4e0e-b54e-81f863aba707",
}
