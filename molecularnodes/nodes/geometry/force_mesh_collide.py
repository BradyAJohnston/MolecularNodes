# Node-group asset "Force Mesh Collide" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    FloatSocket,
    GeometrySocket,
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputFloat, InputGeometry, InputMenu, InputVector


class ForceMeshCollide(AssetGeometryGroup):
    """
    Force Mesh Collide

    Parameters
    ----------
    add : InputVector
        Add
    geometry : InputGeometry
        Geometry
    geometry_bounds : InputMenu | Literal["Original", "Convex Hull"]
        Geometry Bounds
    collision_distance : InputFloat
        Collision Distance

    Inputs
    ------
    i.add : VectorSocket
        Add
    i.geometry : GeometrySocket
        Geometry
    i.geometry_bounds : MenuSocket
        Geometry Bounds
    i.collision_distance : FloatSocket
        Collision Distance

    Outputs
    -------
    o.force : VectorSocket
        Force
    """

    _name = "Force Mesh Collide"
    _asset_name = "Force Mesh Collide"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "INPUT"

    class _Inputs(SocketAccessor):
        add: VectorSocket
        """Add"""
        geometry: GeometrySocket
        """Geometry"""
        geometry_bounds: MenuSocket
        """Geometry Bounds"""
        collision_distance: FloatSocket
        """Collision Distance"""

    class _Outputs(SocketAccessor):
        force: VectorSocket
        """Force"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        add: InputVector = None,
        geometry: InputGeometry = None,
        geometry_bounds: InputMenu | Literal["Original", "Convex Hull"] = "Original",
        collision_distance: InputFloat = 0.1,
    ):
        super().__init__(
            **{
                "Add": add,
                "Geometry": geometry,
                "Geometry Bounds": geometry_bounds,
                "Collision Distance": collision_distance,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        add = tree.inputs.vector(
            "Add",
            (0.0, 0.0, 0.0),
            min_value=-10_000.0,
            max_value=10_000.0,
            hide_value=True,
        )
        geometry = tree.inputs.geometry("Geometry")
        geometry_bounds = tree.inputs.menu("Geometry Bounds", optional_label=True)
        collision_distance = tree.inputs.float(
            "Collision Distance", 0.1, min_value=-10_000.0, max_value=10_000.0
        )
        force = tree.outputs.vector("Force")

        menu_switch = g.MenuSwitch.geometry(
            geometry_bounds,
            {"Original": geometry, "Convex Hull": g.ConvexHull(geometry=geometry)},
        )
        geometry_proximity = g.GeometryProximity(target=menu_switch)
        map_range = geometry_proximity.o.distance.map_range(
            from_max=collision_distance, to_min=1.0, to_max=0.0
        )
        vector_math = g.Position().o.position - geometry_proximity.o.position
        vector_math_1 = vector_math.normalize()
        vector_math_2 = g.Raycast.vector(
            menu_switch, (0.0, 0.0, 0.0), ray_direction=vector_math * -1.0
        ).o.hit_normal.dot(vector_math_1)
        (
            add + vector_math_1 * (vector_math_2 > -0.1).switch.float(-1.0, map_range)
            >> force
        )

        geometry_bounds.default_value = "Original"


ASSET = ForceMeshCollide

ASSET_METADATA = {
    "catalog_id": "c2c958af-5095-4fc2-884d-709bba965fc4",
}
