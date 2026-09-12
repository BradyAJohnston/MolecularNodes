# Node-group asset "Plexus" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    GeometrySocket,
    PackageLibrary,
    SocketAccessor,
)
from nodebpy.types import InputBoolean, InputFloat, InputGeometry
from .sample_position import SamplePosition


class Plexus(AssetGeometryGroup):
    """
    Plexus

    Parameters
    ----------
    points : InputGeometry
        Points
    sort : InputBoolean
        Sort the resulting points so their `Index` matches the input geometry
    distance : InputFloat
        Distance
    radius : InputFloat
        Radius

    Inputs
    ------
    i.points : GeometrySocket
        Points
    i.sort : BooleanSocket
        Sort the resulting points so their `Index` matches the input geometry
    i.distance : FloatSocket
        Distance
    i.radius : FloatSocket
        Radius

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "Plexus"
    _asset_name = "Plexus"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        points: GeometrySocket
        """Points"""
        sort: BooleanSocket
        """Sort the resulting points so their `Index` matches the input geometry"""
        distance: FloatSocket
        """Distance"""
        radius: FloatSocket
        """Radius"""

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
        points: InputGeometry = None,
        sort: InputBoolean = True,
        distance: InputFloat = 0.5,
        radius: InputFloat = 1.0,
    ):
        super().__init__(
            **{"Points": points, "Sort": sort, "Distance": distance, "Radius": radius}
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        points = tree.inputs.geometry("Points")
        sort = tree.inputs.boolean(
            "Sort",
            True,
            description="Sort the resulting points so their `Index` matches the input geometry",
        )
        distance = tree.inputs.float(
            "Distance", 0.5, min_value=-10_000.0, max_value=10_000.0, subtype="DISTANCE"
        )
        radius = tree.inputs.float("Radius", 1.0)
        geometry = tree.outputs.geometry("Geometry")

        with g.Frame("Create a clean set of points for instancing on"):
            index = g.Index()
            sample_index = g.SampleIndex(
                geometry=points,
                value=g.Position(),
                index=index,
                data_type="FLOAT_VECTOR",
            )
            math_1 = g.SampleIndex(
                geometry=points, value=distance, index=index
            ).o.value * g.SampleIndex(geometry=points, value=radius, index=index)
            points_1 = g.Points(
                count=g.DomainSize(geometry=points).o.point_count,
                position=sample_index,
                radius=math_1,
            )
        with g.Frame("Create Distance Probe"):
            ico_sphere = g.IcoSphere()
            sample_index_1 = g.SampleIndex(
                geometry=ico_sphere,
                value=g.Position().o.position * -1.0,
                index=g.Index(),
                data_type="FLOAT_VECTOR",
            )
            merge_by_distance = (
                g.InstanceOnPoints(
                    points=ico_sphere,
                    instance=g.MeshLine(count=2),
                    rotation=g.AxesToRotation(primary_axis=sample_index_1),
                )
                >> g.RealizeInstances(realize_to_point_domain=True)
                >> g.MergeByDistance(distance=0.001)
            )
        with g.Frame("Apply the distance probe"):
            realize_instances = g.RealizeInstances(
                geometry=g.InstanceOnPoints(
                    points=points_1, instance=merge_by_distance, scale=g.Radius()
                ),
                realize_to_point_domain=True,
            )
            capture = g.CaptureAttribute.point(geometry=realize_instances)
            index_1 = capture.items.integer("Index", g.SampleNearest.point(points_1))
            merge_points = (
                capture.o.geometry
                >> g.SetPosition(
                    position=SamplePosition(geometry=points_1, index=index_1.output)
                )
                >> g.MergePoints(merge_id=index_1.output)
            )
        with g.Frame("Potentially sort the results into same order"):
            (
                sort.switch.geometry(
                    merge_points,
                    g.SortElements.point(
                        merge_points, sort_weight=g.SampleNearest.point(points)
                    ),
                )
                >> geometry
            )


ASSET = Plexus

ASSET_METADATA = {
    "catalog_id": "a1e4128a-131f-4e0e-b54e-81f863aba707",
}
