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
    FloatSocket,
    GeometrySocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputClosure, InputFloat, InputGeometry
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
    positive_space : InputBoolean
        Also shift the planar geometry so that all of it has positive coordinates, starting at Padding. Grids built inside the closure then only need OpenVDB tree nodes on one side of each axis, which is faster
    padding : InputFloat
        Distance from the origin to the planar geometry's bounding box when Positive Space is enabled. Should cover anything the closure grows beyond the geometry, such as surface offsets

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry to transform
    i.selection : BooleanSocket
        The parts of the geometry that contribute to the planar calculation
    i.closure : ClosureSocket
        Closure
    i.positive_space : BooleanSocket
        Also shift the planar geometry so that all of it has positive coordinates, starting at Padding. Grids built inside the closure then only need OpenVDB tree nodes on one side of each axis, which is faster
    i.padding : FloatSocket
        Distance from the origin to the planar geometry's bounding box when Positive Space is enabled. Should cover anything the closure grows beyond the geometry, such as surface offsets

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
        positive_space: BooleanSocket
        """Also shift the planar geometry so that all of it has positive coordinates, starting at Padding. Grids built inside the closure then only need OpenVDB tree nodes on one side of each axis, which is faster"""
        padding: FloatSocket
        """Distance from the origin to the planar geometry's bounding box when Positive Space is enabled. Should cover anything the closure grows beyond the geometry, such as surface offsets"""

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
        positive_space: InputBoolean = False,
        padding: InputFloat = 1.0,
    ):
        super().__init__(
            **{
                "Geometry": geometry,
                "Selection": selection,
                "Closure": closure,
                "Positive Space": positive_space,
                "Padding": padding,
            }
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
        positive_space = tree.inputs.boolean(
            "Positive Space",
            False,
            description="Also shift the planar geometry so that all of it has positive coordinates, starting at Padding. Grids built inside the closure then only need OpenVDB tree nodes on one side of each axis, which is faster",
        )
        padding = tree.inputs.float(
            "Padding",
            1.0,
            description="Distance from the origin to the planar geometry's bounding box when Positive Space is enabled. Should cover anything the closure grows beyond the geometry, such as surface offsets",
            min_value=0.0,
            subtype="DISTANCE",
        )
        geometry_1 = tree.outputs.geometry("Geometry")

        geometry_to_planar = GeometryToPlanar(geometry=geometry, selection=selection)
        with g.Frame("Shift into positive space"):
            switch = positive_space.switch.vector(
                (0.0, 0.0, 0.0),
                g.BoundingBox(geometry=geometry_to_planar).o.min * -1.0 + padding,
            )
            transform_geometry = g.TransformGeometry(
                geometry=geometry_to_planar, translation=switch
            )
        evaluate_closure = g.EvaluateClosure(closure)
        evaluate_closure.inputs.geometry(transform_geometry, "Geometry")
        geometry_2 = evaluate_closure.outputs.geometry("Geometry")
        (
            geometry_2.output
            >> g.TransformGeometry(translation=switch * -1.0)
            >> g.TransformGeometry(
                transform=geometry_to_planar.o.transform.invert(), mode="Matrix"
            )
            >> geometry_1
        )


ASSET = EvaluateWhilePlanar

ASSET_METADATA = {
    "catalog_id": "a1e4128a-131f-4e0e-b54e-81f863aba707",
}
