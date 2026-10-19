# Node-group asset "Best View" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
import math
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    FloatSocket,
    GeometrySocket,
    IntegerSocket,
    MatrixSocket,
    PackageLibrary,
    RotationSocket,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import (
    InputBoolean,
    InputFloat,
    InputGeometry,
    InputInteger,
    InputVector,
)
from ._shared.score_view_directions import ScoreViewDirections
from .angstrom_to_world import AngstromToWorld


class BestView(AssetGeometryGroup):
    """
    The direction from which the most selected points have a clear line of sight past the occluding meshes, as a camera transform

    Parameters
    ----------
    target : InputGeometry
        Geometry whose points should be seen, such as the atoms
    selection : InputBoolean
        The points of Target that should be seen, such as a ligand
    occluders : InputGeometry
        Meshes that can block the view, such as the output of a style. Point clouds and instances do not block rays
    up : InputVector
        Direction kept upright in the view
    stickiness : InputFloat
        Share of the selected points a turn of 180° away from Previous Direction must gain to be worth it. Shorter turns cost proportionally less
    previous_direction : InputVector
        Direction currently looked along, which Stickiness favours and which is always scored. Zero for none
    directions : InputInteger
        Directions spread evenly over the sphere in the first pass
    refine : InputInteger
        Rings of refinement around the best direction, starting 10° away and shrinking by 0.4 each time
    offset : InputFloat
        How far in Å towards the camera each ray starts, to clear the point's own geometry

    Inputs
    ------
    i.target : GeometrySocket
        Geometry whose points should be seen, such as the atoms
    i.selection : BooleanSocket
        The points of Target that should be seen, such as a ligand
    i.occluders : GeometrySocket
        Meshes that can block the view, such as the output of a style. Point clouds and instances do not block rays
    i.up : VectorSocket
        Direction kept upright in the view
    i.stickiness : FloatSocket
        Share of the selected points a turn of 180° away from Previous Direction must gain to be worth it. Shorter turns cost proportionally less
    i.previous_direction : VectorSocket
        Direction currently looked along, which Stickiness favours and which is always scored. Zero for none
    i.directions : IntegerSocket
        Directions spread evenly over the sphere in the first pass
    i.refine : IntegerSocket
        Rings of refinement around the best direction, starting 10° away and shrinking by 0.4 each time
    i.offset : FloatSocket
        How far in Å towards the camera each ray starts, to clear the point's own geometry

    Outputs
    -------
    o.transform : MatrixSocket
        Camera transform: located at the centre of the selection and looking down its local -Z along Direction
    o.rotation : RotationSocket
        Camera rotation looking along Direction
    o.direction : VectorSocket
        Direction the camera looks along
    o.centre : VectorSocket
        Centre of the selected points
    o.visible : IntegerSocket
        Selected points with a clear line of sight
    o.count : IntegerSocket
        Selected points
    """

    _name = "Best View"
    _asset_name = "Best View"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"
    _tree_properties = {
        "description": "The direction from which the most selected points have a clear line of sight past the occluding meshes, as a camera transform"
    }

    class _Inputs(SocketAccessor):
        target: GeometrySocket
        """Geometry whose points should be seen, such as the atoms"""
        selection: BooleanSocket
        """The points of Target that should be seen, such as a ligand"""
        occluders: GeometrySocket
        """Meshes that can block the view, such as the output of a style. Point clouds and instances do not block rays"""
        up: VectorSocket
        """Direction kept upright in the view"""
        stickiness: FloatSocket
        """Share of the selected points a turn of 180° away from Previous Direction must gain to be worth it. Shorter turns cost proportionally less"""
        previous_direction: VectorSocket
        """Direction currently looked along, which Stickiness favours and which is always scored. Zero for none"""
        directions: IntegerSocket
        """Directions spread evenly over the sphere in the first pass"""
        refine: IntegerSocket
        """Rings of refinement around the best direction, starting 10° away and shrinking by 0.4 each time"""
        offset: FloatSocket
        """How far in Å towards the camera each ray starts, to clear the point's own geometry"""

    class _Outputs(SocketAccessor):
        transform: MatrixSocket
        """Camera transform: located at the centre of the selection and looking down its local -Z along Direction"""
        rotation: RotationSocket
        """Camera rotation looking along Direction"""
        direction: VectorSocket
        """Direction the camera looks along"""
        centre: VectorSocket
        """Centre of the selected points"""
        visible: IntegerSocket
        """Selected points with a clear line of sight"""
        count: IntegerSocket
        """Selected points"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        target: InputGeometry = None,
        selection: InputBoolean = True,
        occluders: InputGeometry = None,
        up: InputVector = None,
        stickiness: InputFloat = 0.0,
        previous_direction: InputVector = None,
        directions: InputInteger = 256,
        refine: InputInteger = 3,
        offset: InputFloat = 1.7,
    ):
        super().__init__(
            **{
                "Target": target,
                "Selection": selection,
                "Occluders": occluders,
                "Up": up,
                "Stickiness": stickiness,
                "Previous Direction": previous_direction,
                "Directions": directions,
                "Refine": refine,
                "Offset": offset,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        target = tree.inputs.geometry(
            "Target",
            description="Geometry whose points should be seen, such as the atoms",
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="The points of Target that should be seen, such as a ligand",
            hide_value=True,
        )
        occluders = tree.inputs.geometry(
            "Occluders",
            description="Meshes that can block the view, such as the output of a style. Point clouds and instances do not block rays",
        )
        up = tree.inputs.vector(
            "Up",
            (0.0, 0.0, 1.0),
            description="Direction kept upright in the view",
            subtype="XYZ",
        )
        stickiness = tree.inputs.float(
            "Stickiness",
            0.0,
            description="Share of the selected points a turn of 180° away from Previous Direction must gain to be worth it. Shorter turns cost proportionally less",
            min_value=0.0,
            max_value=1.0,
            subtype="FACTOR",
        )
        previous_direction = tree.inputs.vector(
            "Previous Direction",
            (0.0, 0.0, 0.0),
            description="Direction currently looked along, which Stickiness favours and which is always scored. Zero for none",
        )
        with tree.inputs.panel("Search", default_closed=True):
            directions = tree.inputs.integer(
                "Directions",
                256,
                description="Directions spread evenly over the sphere in the first pass",
                min_value=8,
                max_value=4096,
            )
            refine = tree.inputs.integer(
                "Refine",
                3,
                description="Rings of refinement around the best direction, starting 10° away and shrinking by 0.4 each time",
                min_value=0,
                max_value=8,
            )
            offset = tree.inputs.float(
                "Offset",
                1.7,
                description="How far in Å towards the camera each ray starts, to clear the point's own geometry",
                min_value=0.0,
            )
        transform = tree.outputs.matrix(
            "Transform",
            description="Camera transform: located at the centre of the selection and looking down its local -Z along Direction",
        )
        rotation = tree.outputs.rotation(
            "Rotation", description="Camera rotation looking along Direction"
        )
        direction = tree.outputs.vector(
            "Direction", description="Direction the camera looks along"
        )
        centre = tree.outputs.vector(
            "Centre", description="Centre of the selected points"
        )
        visible = tree.outputs.integer(
            "Visible", description="Selected points with a clear line of sight"
        )
        count = tree.outputs.integer("Count", description="Selected points")

        with g.Frame("Selected points only"):
            angstrom_to_world = AngstromToWorld(angstrom=offset)
            separate_geometry = g.SeparateGeometry.point(target, selection)
            separate_components = g.SeparateComponents(separate_geometry.o.selection)
            join_geometry = g.JoinGeometry(
                (
                    g.MeshToPoints(separate_geometry.o.selection),
                    separate_components.o.point_cloud,
                    g.InstancesToPoints(separate_components.o.instances),
                )
            )
            domain_size = g.DomainSize(geometry=join_geometry, component="POINTCLOUD")
            points = g.Points(
                count=domain_size.o.point_count,
                position=g.SampleIndex.point.vector(
                    join_geometry, g.Position(), g.Index()
                ),
            )
        with g.Frame("Fibonacci sphere of directions"):
            integer_math = directions.max(8)
            math_1 = g.Index().o.index + 0.5
            math_2 = 1.0 - math_1 * 2.0 / integer_math
            math_3 = (1.0 - math_2 * math_2).max(0.0).sqrt()
            math_4 = math_1 * 2.399963
            _string = g.String(
                "A Previous Direction is added as one more candidate, so staying put is always scored exactly rather than at the nearest point of the sphere."
            )
            switch = (g.Index() < integer_math).switch.vector(
                previous_direction.normalize(),
                g.CombineXYZ(
                    x=math_3 * math_4.cos(), y=math_3 * math_4.sin(), z=math_2
                ),
            )
            points_1 = g.Points(
                count=integer_math
                + (previous_direction.length() > 0.5).switch.integer(true=1),
                position=switch,
            )
            score_view_directions = ScoreViewDirections(
                directions=points_1,
                target=points,
                occluders=occluders,
                offset=angstrom_to_world,
                previous=previous_direction,
                stickiness=stickiness,
            )
        with g.Frame("Refine on shrinking rings"):
            math_5 = g.Index().o.index * (math.pi / 8)
            _string_1 = g.String(
                "Scores 16 directions on a ring around the current best and moves to the best of them only if it beats the current score."
            )
            repeat_zone = g.RepeatZone(refine)
            direction_1 = repeat_zone.items.vector(
                score_view_directions.o.direction, "Direction"
            )
            score = repeat_zone.items.float(score_view_directions.o.score, "Score")
            visible_1 = repeat_zone.items.integer(
                score_view_directions.o.visible, "Visible"
            )
            math_6 = 0.4**repeat_zone.iteration * 0.1745329
            vector_math = (
                (abs(direction_1.current.x) < 0.9)
                .switch.vector((0.0, 1.0, 0.0), (1.0, 0.0, 0.0))
                .cross(direction_1.current)
            )
            vector_math_1 = vector_math.normalize()
            vector_math_2 = (
                vector_math_1 * math_5.cos()
                + direction_1.current.cross(vector_math_1) * math_5.sin()
            ) * math_6.sin()
            points_2 = g.Points(
                position=(
                    direction_1.current * math_6.cos() + vector_math_2
                ).normalize(),
                count=16,
            )
            score_view_directions_1 = ScoreViewDirections(
                directions=points_2,
                target=points,
                occluders=occluders,
                offset=angstrom_to_world,
                previous=previous_direction,
                stickiness=stickiness,
            )
            compare = score_view_directions_1.o.score > score.current
            (
                compare.switch.vector(
                    direction_1.current, score_view_directions_1.o.direction
                )
                >> direction_1.next
            )
            (
                compare.switch.float(score.current, score_view_directions_1.o.score)
                >> score.next
            )
            (
                compare.switch.integer(
                    visible_1.current, score_view_directions_1.o.visible
                )
                >> visible_1.next
            )
        with g.Frame("Camera looking along the direction"):
            sample_index = g.SampleIndex.point.vector(
                points, g.Position().o.position.point.mean()
            )
            vector_math_3 = up.normalize()
            switch_1 = (
                abs(direction_1.result.dot(vector_math_3)) > 0.95
            ).switch.vector(vector_math_3, (0.0, 1.0, 0.0))
            axes_to_rotation = g.AxesToRotation(
                primary_axis=direction_1.result * -1.0,
                secondary_axis=switch_1,
                secondary="Y",
            )
            combine_transform = g.CombineTransform(
                translation=sample_index, rotation=axes_to_rotation
            )

        combine_transform >> transform
        axes_to_rotation >> rotation
        direction_1.result >> direction
        sample_index >> centre
        visible_1.result >> visible
        domain_size >> count


ASSET = BestView

ASSET_METADATA = {
    "description": "The direction from which the most selected points have a clear line of sight past the occluding meshes, as a camera transform",
    "catalog_id": "a1e4128a-131f-4e0e-b54e-81f863aba707",
}
