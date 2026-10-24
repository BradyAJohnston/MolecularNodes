# Node group "Score View Directions" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (nodebpy build).
# Shared by several assets, which import it; not an asset itself.
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    CustomGeometryGroup,
    FloatSocket,
    GeometrySocket,
    IntegerSocket,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import InputFloat, InputGeometry, InputVector


class ScoreViewDirections(CustomGeometryGroup):
    """
    Score View Directions

    Parameters
    ----------
    directions : InputGeometry
        Points whose positions are the unit directions the camera could look along
    target : InputGeometry
        Points that should be seen from the camera
    occluders : InputGeometry
        Meshes that can block the line of sight
    offset : InputFloat
        How far towards the camera each ray starts, to clear the point's own geometry
    previous : InputVector
        Direction currently looked along. Zero for none
    stickiness : InputFloat
        Share of the target points a turn of 180° away from Previous must gain to be worth it. Shorter turns cost proportionally less

    Inputs
    ------
    i.directions : GeometrySocket
        Points whose positions are the unit directions the camera could look along
    i.target : GeometrySocket
        Points that should be seen from the camera
    i.occluders : GeometrySocket
        Meshes that can block the line of sight
    i.offset : FloatSocket
        How far towards the camera each ray starts, to clear the point's own geometry
    i.previous : VectorSocket
        Direction currently looked along. Zero for none
    i.stickiness : FloatSocket
        Share of the target points a turn of 180° away from Previous must gain to be worth it. Shorter turns cost proportionally less

    Outputs
    -------
    o.direction : VectorSocket
        The best scoring direction
    o.score : FloatSocket
        Visible points of the best direction, plus a fraction below one that favours a broader projection
    o.visible : IntegerSocket
        Target points with a clear line of sight
    """

    _name = "Score View Directions"
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        directions: GeometrySocket
        """Points whose positions are the unit directions the camera could look along"""
        target: GeometrySocket
        """Points that should be seen from the camera"""
        occluders: GeometrySocket
        """Meshes that can block the line of sight"""
        offset: FloatSocket
        """How far towards the camera each ray starts, to clear the point's own geometry"""
        previous: VectorSocket
        """Direction currently looked along. Zero for none"""
        stickiness: FloatSocket
        """Share of the target points a turn of 180° away from Previous must gain to be worth it. Shorter turns cost proportionally less"""

    class _Outputs(SocketAccessor):
        direction: VectorSocket
        """The best scoring direction"""
        score: FloatSocket
        """Visible points of the best direction, plus a fraction below one that favours a broader projection"""
        visible: IntegerSocket
        """Target points with a clear line of sight"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        directions: InputGeometry = None,
        target: InputGeometry = None,
        occluders: InputGeometry = None,
        offset: InputFloat = 0.17,
        previous: InputVector = None,
        stickiness: InputFloat = 0.0,
    ):
        super().__init__(
            **{
                "Directions": directions,
                "Target": target,
                "Occluders": occluders,
                "Offset": offset,
                "Previous": previous,
                "Stickiness": stickiness,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        directions = tree.inputs.geometry(
            "Directions",
            description="Points whose positions are the unit directions the camera could look along",
        )
        target = tree.inputs.geometry(
            "Target", description="Points that should be seen from the camera"
        )
        occluders = tree.inputs.geometry(
            "Occluders", description="Meshes that can block the line of sight"
        )
        offset = tree.inputs.float(
            "Offset",
            0.17,
            description="How far towards the camera each ray starts, to clear the point's own geometry",
            min_value=0.0,
            subtype="DISTANCE",
        )
        previous = tree.inputs.vector(
            "Previous",
            (0.0, 0.0, 0.0),
            description="Direction currently looked along. Zero for none",
        )
        stickiness = tree.inputs.float(
            "Stickiness",
            0.0,
            description="Share of the target points a turn of 180° away from Previous must gain to be worth it. Shorter turns cost proportionally less",
            min_value=0.0,
        )
        direction = tree.outputs.vector(
            "Direction", description="The best scoring direction"
        )
        score = tree.outputs.float(
            "Score",
            description="Visible points of the best direction, plus a fraction below one that favours a broader projection",
        )
        visible = tree.outputs.integer(
            "Visible", description="Target points with a clear line of sight"
        )

        domain_size = g.DomainSize(geometry=target, component="POINTCLOUD")
        with g.Frame("Every direction paired with every target point"):
            _string = g.String(
                "Each direction gets its own copy of the target points, so one Raycast scores all directions at once. The candidate index groups the copies back together."
            )
            capture = g.CaptureAttribute.point(geometry=directions)
            forward = capture.items.vector(
                g.Position().o.position.normalize(), "Forward"
            )
            candidate = capture.items.integer(g.Index(), "Candidate")
            realize_instances = (
                g.SetPosition(
                    geometry=capture.o.geometry, position=g.Vector((0.0, 0.0, 0.0))
                )
                >> g.InstanceOnPoints(instance=target)
                >> g.RealizeInstances()
            )
        with g.Frame("Rays from each point back towards the camera"):
            position = g.Position()
            vector_math = forward.output * -1.0
            raycast = g.Raycast.float(
                occluders,
                source_position=position.o.position + vector_math * offset,
                ray_direction=vector_math,
                ray_length=1000.0,
            )
            accumulate_field = raycast.o.is_hit.switch.integer(1).point.total(
                candidate.output
            )
        with g.Frame("Tie-break on projected extent"):
            _string_1 = g.String(
                "Adds spread / (spread + 1) / 2 for the projected bounding box area, so among equally clear views the one showing the target broad side on wins without ever outweighing a visible point."
            )
            vector_math_1 = (
                (abs(forward.output.x) < 0.9)
                .switch.vector((0.0, 1.0, 0.0), (1.0, 0.0, 0.0))
                .cross(forward.output)
            )
            vector_math_2 = vector_math_1.normalize()
            vector_math_3 = position.o.position.dot(vector_math_2)
            vector_math_4 = position.o.position.dot(forward.output.cross(vector_math_2))
            math_1 = (
                vector_math_3.point.max(candidate.output)
                - vector_math_3.point.min(candidate.output)
            ) * (
                vector_math_4.point.max(candidate.output)
                - vector_math_4.point.min(candidate.output)
            )
            capture_1 = g.CaptureAttribute.point(geometry=realize_instances)
            score_1 = capture_1.items.float(
                accumulate_field + math_1 / (math_1 + 1.0) * 0.5, "Score"
            )
            visible_1 = capture_1.items.integer(accumulate_field, "Visible")
        with g.Frame("Pick the best direction"):
            integer_math = g.Index().o.index * domain_size.o.point_count
            _string_2 = g.String(
                "Stickiness subtracts its share of the points scaled by (1 - cos) / 2 of the angle from Previous, so staying put is free and turning around costs the whole share. With Previous at zero every direction pays the same and the choice is unchanged."
            )
            math_2 = (
                g.SampleIndex.point.float(
                    capture_1.o.geometry, score_1.output, integer_math
                ).o.value
                - stickiness
                * domain_size.o.point_count
                * (1.0 - forward.output.dot(previous))
                * 0.5
            )
            capture_2 = g.CaptureAttribute.point(geometry=capture.o.geometry)
            score_2 = capture_2.items.float(math_2, "Score")
            visible_2 = capture_2.items.integer(
                g.SampleIndex.point.integer(
                    capture_1.o.geometry, visible_1.output, integer_math
                ),
                "Visible",
            )
            field_min_max = (
                (score_2.output >= score_2.output.point.max())
                .switch.integer(1073741824, g.Index())
                .point.min()
            )
            sample_index = g.SampleIndex.point.integer(
                capture_2.o.geometry, field_min_max
            )
            sample_index_1 = g.SampleIndex.point.vector(
                capture_2.o.geometry, forward.output, sample_index
            )
            sample_index_2 = g.SampleIndex.point.float(
                capture_2.o.geometry, score_2.output, sample_index
            )
            sample_index_3 = g.SampleIndex.point.integer(
                capture_2.o.geometry, visible_2.output, sample_index
            )

        sample_index_1 >> direction
        sample_index_2 >> score
        sample_index_3 >> visible
