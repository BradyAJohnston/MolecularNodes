# Node-group asset "Evluate While Planar" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
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
from .geoemtry_to_planar import GeoemtryToPlanar


class EvluateWhilePlanar(AssetGeometryGroup):
    """
    Evluate While Planar

    Parameters
    ----------
    geometry : InputGeometry
        Geometry to transform
    selection : InputBoolean
        The parts of the geometry that contibute to the planar calculation
    closure : InputClosure
        Closure

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry to transform
    i.selection : BooleanSocket
        The parts of the geometry that contibute to the planar calculation
    i.closure : ClosureSocket
        Closure

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "Evluate While Planar"
    _asset_name = "Evluate While Planar"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry to transform"""
        selection: BooleanSocket
        """The parts of the geometry that contibute to the planar calculation"""
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
            description="The parts of the geometry that contibute to the planar calculation",
            hide_value=True,
        )
        closure = tree.inputs.closure("Closure")
        geometry_1 = tree.outputs.geometry("Geometry")

        group = GeoemtryToPlanar(geometry=geometry, selection=selection)
        evaluate_closure = g.EvaluateClosure(closure)
        evaluate_closure.inputs.geometry("Geometry", group.o.geometry)
        geometry_2 = evaluate_closure.outputs.geometry("Geometry")
        (
            geometry_2
            >> g.TransformGeometry(transform=group.o.transform.invert(), mode="Matrix")
            >> geometry_1
        )


ASSET = EvluateWhilePlanar

ASSET_METADATA = {
    "catalog_id": "a1e4128a-131f-4e0e-b54e-81f863aba707",
}
