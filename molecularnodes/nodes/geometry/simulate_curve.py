# Node-group asset "Simulate Curve" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING
from bpy.types import GeometryNodeTree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder import (
    AssetGeometryGroup,
    BooleanSocket,
    CustomGeometryGroup,
    FloatSocket,
    GeometrySocket,
    IntegerSocket,
    PackageLibrary,
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
from ._shared.constraint_distance import ConstraintDistance
from ._shared.inverse_mass import InverseMass
from ._shared.xpbd_finalise import XPBDFinalise
from ._shared.xpbd_init import XPBDInit
from ._shared.xpbd_solve_hook import XPBDSolveHook
from .mass import Mass
from .offset_float import OffsetFloat
from .offset_vector import OffsetVector


class XPBDSolveCurve(CustomGeometryGroup):
    _name = "XPBD Solve Curve"
    _color_tag = "GEOMETRY"

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry("Geometry")
        selection = tree.inputs.boolean("Selection", True, hide_value=True)
        points = tree.inputs.integer("Points", 2, min_value=0, max_value=5)
        straightness = tree.inputs.float(
            "Straightness", 0.9, min_value=0.0, max_value=1.0, subtype="FACTOR"
        )
        length = tree.inputs.float("Length", 0.0)
        alpha = tree.inputs.float("alpha", 0.0, min_value=0.0, max_value=10_000.0)
        deltat = tree.inputs.float(
            "deltaT",
            0.0,
            min_value=-10_000.0,
            max_value=10_000.0,
            structure_type="SINGLE",
            force_non_field=True,
        )
        geometry_1 = tree.outputs.geometry("Geometry")

        group = InverseMass()
        separate_components = g.SeparateComponents(geometry=geometry)
        repeat_zone = g.RepeatZone(points)
        geometry_2 = repeat_zone.items.geometry("Geometry", separate_components.o.curve)
        integer_math = repeat_zone.iteration + 1
        math_1 = (
            integer_math
            * g.Mix(
                factor_float=straightness, a_float=0.5, b_float=1.0, clamp_factor=True
            ).o.result_float
        )
        repeat_zone_1 = g.RepeatZone(2)
        geometry_3 = repeat_zone_1.items.geometry("Geometry", geometry_2.current)
        index_switch = g.IndexSwitch.boolean(
            repeat_zone_1.iteration,
            (
                g.EndpointSelection(end_size=integer_math, start_size=0),
                g.EndpointSelection(start_size=integer_math, end_size=0),
            ),
        )
        index_switch_1 = g.IndexSwitch.integer(
            repeat_zone_1.iteration, (integer_math, -integer_math)
        )
        group_1 = ConstraintDistance(
            target=OffsetVector(offset=index_switch_1),
            distance=(repeat_zone.iteration > 0).switch.float(integer_math, math_1)
            * length,
            w1=group.o.w,
            w2=OffsetFloat(value=group.o.w, offset=index_switch_1),
            alpha=alpha,
            deltat=deltat,
        )
        set_position = g.SetPosition(
            geometry=geometry_3.current,
            selection=g.BooleanMath.subtract(selection, index_switch),
            offset=group_1.o.correction,
        )
        set_position >> geometry_3.next
        geometry_3.result >> geometry_2.next
        join_geometry = g.JoinGeometry(
            geometry=(
                separate_components.o.mesh,
                geometry_2.result,
                separate_components.o.grease_pencil,
                separate_components.o.point_cloud,
                separate_components.o.volume,
                separate_components.o.instances,
            )
        )

        join_geometry >> geometry_1


class XPBDSolvePointsForCurve(CustomGeometryGroup):
    _name = "XPBD Solve Points for Curve"
    _color_tag = "GEOMETRY"

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry("Geometry")
        selection = tree.inputs.boolean("Selection", True, hide_value=True)
        radius = tree.inputs.float("Radius", 0.0)
        alpha = tree.inputs.float("alpha", 0.0, min_value=-10_000.0, max_value=10_000.0)
        deltat = tree.inputs.float(
            "deltaT",
            0.0,
            min_value=-10_000.0,
            max_value=10_000.0,
            structure_type="SINGLE",
            force_non_field=True,
        )
        geometry_1 = tree.outputs.geometry("Geometry")

        vector_math = (
            g.Position().o.position
            + g.RandomValue.vector((-1.0, -1.0, -1.0), (1.0, 1.0, 1.0)).o.value * 0.001
        )
        capture = g.CaptureAttribute.point(geometry=geometry)
        vector = capture.items.vector("Vector", vector_math)
        capture_1 = g.CaptureAttribute.point(geometry=capture.o.geometry)
        index = capture_1.items.integer(
            "Index", g.IndexOfNearest(position=vector.output).o.index
        )
        math_1 = radius.point.at(index.output) + radius
        group = ConstraintDistance(
            target=g.Position().o.position.point.at(index.output),
            distance=math_1,
            w1=1.0,
            w2=1.0,
            alpha=alpha,
            deltat=deltat,
        )
        (
            capture_1.o.geometry
            >> g.SetPosition(
                selection=selection & (math_1 > group.o.value),
                offset=group.o.correction,
            )
            >> geometry_1
        )


class SimulateCurve(AssetGeometryGroup):
    """
    Simulate Curve

    Parameters
    ----------
    geometry : InputGeometry
        Geometry containing curve to simulate
    selection : InputBoolean
        Which points to simulate on the curve
    substeps : InputInteger
        Number of substeps to simulate for each frame of the animation (higher is slower but more accurate)
    force : InputVector
        Additional forces to be applied to the simulation
    drag : InputFloat
        Drag force for how quickly velocity should dissipate in the simulation. 0 is no drag and higher values dissipate quicker.
    point_radius : InputFloat
        Radius for each point which is used for calculating collisions between points
    point_alpha : InputFloat
        Point alpha
    curve_points : InputInteger
        Number of points to check in each direction to maintain straightness of curve
    curve_straightness : InputFloat
        Curve Straightness
    curve_segment_length : InputFloat
        Curve Segment Length
    curve_alpha : InputFloat
        Curve alpha
    hook_selection : InputBoolean
        Hook Selection
    hook_target : InputVector
        Hook Target
    hook_decay : InputFloat
        Hook Decay
    pin_selection : InputBoolean
        Pin Selection
    pin_target : InputVector
        Pin Target

    Inputs
    ------
    i.geometry : GeometrySocket
        Geometry containing curve to simulate
    i.selection : BooleanSocket
        Which points to simulate on the curve
    i.substeps : IntegerSocket
        Number of substeps to simulate for each frame of the animation (higher is slower but more accurate)
    i.force : VectorSocket
        Additional forces to be applied to the simulation
    i.drag : FloatSocket
        Drag force for how quickly velocity should dissipate in the simulation. 0 is no drag and higher values dissipate quicker.
    i.point_radius : FloatSocket
        Radius for each point which is used for calculating collisions between points
    i.point_alpha : FloatSocket
        Point alpha
    i.curve_points : IntegerSocket
        Number of points to check in each direction to maintain straightness of curve
    i.curve_straightness : FloatSocket
        Curve Straightness
    i.curve_segment_length : FloatSocket
        Curve Segment Length
    i.curve_alpha : FloatSocket
        Curve alpha
    i.hook_selection : BooleanSocket
        Hook Selection
    i.hook_target : VectorSocket
        Hook Target
    i.hook_decay : FloatSocket
        Hook Decay
    i.pin_selection : BooleanSocket
        Pin Selection
    i.pin_target : VectorSocket
        Pin Target

    Outputs
    -------
    o.geometry : GeometrySocket
        Geometry
    """

    _name = "Simulate Curve"
    _asset_name = "Simulate Curve"
    _library = PackageLibrary(__file__, "../../assets/node_data_file.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        geometry: GeometrySocket
        """Geometry containing curve to simulate"""
        selection: BooleanSocket
        """Which points to simulate on the curve"""
        substeps: IntegerSocket
        """Number of substeps to simulate for each frame of the animation (higher is slower but more accurate)"""
        force: VectorSocket
        """Additional forces to be applied to the simulation"""
        drag: FloatSocket
        """Drag force for how quickly velocity should dissipate in the simulation. 0 is no drag and higher values dissipate quicker."""
        point_radius: FloatSocket
        """Radius for each point which is used for calculating collisions between points"""
        point_alpha: FloatSocket
        """Point alpha"""
        curve_points: IntegerSocket
        """Number of points to check in each direction to maintain straightness of curve"""
        curve_straightness: FloatSocket
        """Curve Straightness"""
        curve_segment_length: FloatSocket
        """Curve Segment Length"""
        curve_alpha: FloatSocket
        """Curve alpha"""
        hook_selection: BooleanSocket
        """Hook Selection"""
        hook_target: VectorSocket
        """Hook Target"""
        hook_decay: FloatSocket
        """Hook Decay"""
        pin_selection: BooleanSocket
        """Pin Selection"""
        pin_target: VectorSocket
        """Pin Target"""

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
        substeps: InputInteger = 5,
        force: InputVector = None,
        drag: InputFloat = 1.0,
        point_radius: InputFloat = 0.04,
        point_alpha: InputFloat = 0.0,
        curve_points: InputInteger = 1,
        curve_straightness: InputFloat = 0.9,
        curve_segment_length: InputFloat = 0.1,
        curve_alpha: InputFloat = 0.0,
        hook_selection: InputBoolean = False,
        hook_target: InputVector = None,
        hook_decay: InputFloat = 0.5,
        pin_selection: InputBoolean = False,
        pin_target: InputVector = None,
    ):
        super().__init__(
            **{
                "Geometry": geometry,
                "Selection": selection,
                "Substeps": substeps,
                "Force": force,
                "Drag": drag,
                "Point Radius": point_radius,
                "Point alpha": point_alpha,
                "Curve Points": curve_points,
                "Curve Straightness": curve_straightness,
                "Curve Segment Length": curve_segment_length,
                "Curve alpha": curve_alpha,
                "Hook Selection": hook_selection,
                "Hook Target": hook_target,
                "Hook Decay": hook_decay,
                "Pin Selection": pin_selection,
                "Pin Target": pin_target,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry(
            "Geometry", description="Geometry containing curve to simulate"
        )
        selection = tree.inputs.boolean(
            "Selection",
            True,
            description="Which points to simulate on the curve",
            hide_value=True,
        )
        substeps = tree.inputs.integer(
            "Substeps",
            5,
            description="Number of substeps to simulate for each frame of the animation (higher is slower but more accurate)",
            min_value=1,
            max_value=100,
        )
        with tree.inputs.panel("Simulation"):
            force = tree.inputs.vector(
                "Force",
                (0.0, 0.0, 0.0),
                description="Additional forces to be applied to the simulation",
                min_value=-10_000.0,
                max_value=10_000.0,
                hide_value=True,
            )
            drag = tree.inputs.float(
                "Drag",
                1.0,
                description="Drag force for how quickly velocity should dissipate in the simulation. 0 is no drag and higher values dissipate quicker.",
                min_value=0.0,
                max_value=10_000.0,
            )
        with tree.inputs.panel("Point"):
            point_radius = tree.inputs.float(
                "Point Radius",
                0.04,
                description="Radius for each point which is used for calculating collisions between points",
                min_value=0.0,
            )
            point_alpha = tree.inputs.float(
                "Point alpha", 0.0, min_value=-10_000.0, max_value=10_000.0
            )
        with tree.inputs.panel("Curve"):
            curve_points = tree.inputs.integer(
                "Curve Points",
                1,
                description="Number of points to check in each direction to maintain straightness of curve",
                min_value=0,
                max_value=5,
            )
            curve_straightness = tree.inputs.float(
                "Curve Straightness",
                0.9,
                min_value=0.0,
                max_value=1.0,
                subtype="FACTOR",
            )
            curve_segment_length = tree.inputs.float(
                "Curve Segment Length", 0.1, min_value=0.0
            )
            curve_alpha = tree.inputs.float(
                "Curve alpha", 0.0, min_value=0.0, max_value=10_000.0
            )
        with tree.inputs.panel("Hook"):
            hook_selection = tree.inputs.boolean(
                "Hook Selection", False, hide_value=True
            )
            hook_target = tree.inputs.vector(
                "Hook Target",
                (0.0, 0.0, 0.0),
                min_value=-10_000.0,
                max_value=10_000.0,
                hide_value=True,
            )
            hook_decay = tree.inputs.float(
                "Hook Decay", 0.5, min_value=-10_000.0, max_value=10_000.0
            )
        with tree.inputs.panel("Pin"):
            pin_selection = tree.inputs.boolean("Pin Selection", False, hide_value=True)
            pin_target = tree.inputs.vector(
                "Pin Target", (0.0, 0.0, 0.0), hide_value=True, default_input="POSITION"
            )
        geometry_1 = tree.outputs.geometry("Geometry")

        boolean_math = g.BooleanMath.subtract(selection, pin_selection)
        store_named_attribute = g.StoreNamedAttribute.point.float(
            g.SeparateComponents(geometry=geometry).o.curve,
            ~g.NamedAttribute.float("mass").o.exists,
            "mass",
            1.0,
        )
        simulation_zone = g.SimulationZone()
        geometry_2 = simulation_zone.items.geometry("Geometry", store_named_attribute)
        math_1 = simulation_zone.delta_time / substeps
        store_named_attribute_1 = g.StoreNamedAttribute.point.float(
            geometry_2.current, name="inverse_mass", value=1.0 / Mass()
        )
        repeat_zone = g.RepeatZone(substeps)
        geometry_3 = repeat_zone.items.geometry("Geometry", store_named_attribute_1)
        group = XPBDInit(
            geometry=geometry_3.current,
            selection=boolean_math,
            force=force,
            drag=drag,
            deltat=math_1,
        )
        group_1 = XPBDSolveCurve(
            Geometry=group,
            Selection=boolean_math,
            Points=curve_points,
            Straightness=curve_straightness,
            Length=curve_segment_length,
            alpha=curve_alpha,
            deltaT=math_1,
        )
        group_2 = XPBDSolvePointsForCurve(
            Geometry=group_1,
            Selection=boolean_math,
            Radius=point_radius,
            alpha=point_alpha,
            deltaT=math_1,
        )
        set_position = XPBDSolveHook(
            geometry=group_2,
            selection=hook_selection,
            target=hook_target,
            decay=hook_decay,
            deltat=math_1,
        ) >> g.SetPosition(selection=pin_selection, position=pin_target)
        (
            XPBDFinalise(geometry=set_position, selection=boolean_math, deltat=math_1)
            >> geometry_3.next
        )
        geometry_3.result >> geometry_2.next

        geometry_2.result >> geometry_1


ASSET = SimulateCurve

ASSET_METADATA = {
    "catalog_id": "c2c958af-5095-4fc2-884d-709bba965fc4",
}
