# Node-group asset 'Lag Geometry' (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    GeometrySocket,
    IntegerSocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputGeometry, InputInteger


class LagGeometry(AssetGeometryGroup):
    """
    Lag Geometry

    Parameters
    ----------
    input : InputGeometry
        Input
    selection : InputBoolean
        The parts of the geometry that go into the first output
    count : InputInteger
        Number of frames to lag.
    realize_all : InputBoolean
        Realize all levels of nested instances for a top-level instances. Overrides the value of the Depth input

    Inputs
    ------
    i.input : GeometrySocket
        Input
    i.selection : BooleanSocket
        The parts of the geometry that go into the first output
    i.count : IntegerSocket
        Number of frames to lag.
    i.realize_all : BooleanSocket
        Realize all levels of nested instances for a top-level instances. Overrides the value of the Depth input

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    o.lag_index : IntegerSocket
        Lag Index
    o.index : IntegerSocket
        Index
    """

    _name = "Lag Geometry"
    _asset_name = "Lag Geometry"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        input: GeometrySocket
        """Input"""
        selection: BooleanSocket
        """The parts of the geometry that go into the first output"""
        count: IntegerSocket
        """Number of frames to lag."""
        realize_all: BooleanSocket
        """Realize all levels of nested instances for a top-level instances. Overrides the value of the Depth input"""

    class _Outputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry"""
        lag_index: IntegerSocket
        """Lag Index"""
        index: IntegerSocket
        """Index"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        input: InputGeometry = None,
        selection: InputBoolean = True,
        count: InputInteger = 5,
        realize_all: InputBoolean = True,
    ):
        super().__init__(
            **{
                "Input": input,
                "Selection": selection,
                "Count": count,
                "Realize All": realize_all,
            }
        )

    def _build_group(self, tree):
        input = tree.inputs.geometry("Input")
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="The parts of the geometry that go into the first output",
            hide_value=True,
        )
        count = tree.inputs.integer("Count", 5, description="Number of frames to lag.")
        realize_all = tree.inputs.boolean(
            "Realize All",
            True,
            description="Realize all levels of nested instances for a top-level instances. Overrides the value of the Depth input",
        )
        geometry = tree.outputs.geometry("Geometry")
        lag_index = tree.outputs.integer("Lag Index")
        index = tree.outputs.integer("Index")

        index_1 = g.Index()
        simulation_zone = g.SimulationZone()
        geometry_1 = simulation_zone.items.geometry("Geometry")
        join_geometry = g.JoinGeometry(
            geometry=(
                g.GeometryToInstance(
                    g.SeparateGeometry.point(input, selection).o.selection
                ),
                geometry_1.current,
            )
        )
        g.SeparateGeometry.instance(join_geometry, g.Index() < count) >> geometry_1.next
        capture = g.CaptureAttribute.point(geometry=geometry_1.result)
        index_2 = capture.items.integer("Index", index_1)
        capture_1 = g.CaptureAttribute.instance(geometry=capture.o.geometry)
        index_3 = capture_1.items.integer("Index", index_1)
        capture_1.o.geometry >> g.RealizeInstances(realize_all=realize_all) >> geometry

        index_3.output >> lag_index
        index_2.output >> index


ASSET = LagGeometry

ASSET_METADATA = {
    "catalog_id": "a1e4128a-131f-4e0e-b54e-81f863aba707",
}
