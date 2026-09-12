# Node-group asset "Simulate on Faces" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
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
from ._shared.xpbd_finalise import XPBDFinalise
from ._shared.xpbd_init import XPBDInit
from ._shared.xpbd_solve_hook import XPBDSolveHook
from ._shared.xpbd_solve_points import XPBDSolvePoints
from .mass import Mass


class XPBDSolveOnFaces(CustomGeometryGroup):
    _name = "XPBD Solve on Faces"
    _color_tag = "GEOMETRY"

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry("Geometry")
        selection = tree.inputs.boolean("Selection", True, hide_value=True)
        faces = tree.inputs.geometry("Faces")
        alpha = tree.inputs.float("alpha", 0.0, min_value=-10_000.0, max_value=10_000.0)
        deltat = tree.inputs.float(
            "deltaT", 0.0, min_value=-10_000.0, max_value=10_000.0
        )
        geometry_1 = tree.outputs.geometry("Geometry")

        geometry_proximity = g.GeometryProximity(target=faces)
        group = ConstraintDistance(
            target=geometry_proximity.o.position,
            distance=0.0,
            alpha=alpha,
            deltat=deltat,
        )
        (
            geometry
            >> g.SetPosition(
                selection=selection & geometry_proximity.o.is_valid,
                offset=group.o.correction,
            )
            >> geometry_1
        )


class SimulateOnFaces(AssetGeometryGroup):
    """
    Simulate on Faces

    Parameters
    ----------
    points : InputGeometry
        Points
    selection : InputBoolean
        Selection
    substeps : InputInteger
        Substeps
    force : InputVector
        Force
    drag : InputFloat
        Drag
    faces : InputGeometry
        Faces
    alpha : InputFloat
        alpha
    particle_radius : InputFloat
        Particle Radius
    particle_alpha : InputFloat
        Particle alpha
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
    i.points : GeometrySocket
        Points
    i.selection : BooleanSocket
        Selection
    i.substeps : IntegerSocket
        Substeps
    i.force : VectorSocket
        Force
    i.drag : FloatSocket
        Drag
    i.faces : GeometrySocket
        Faces
    i.alpha : FloatSocket
        alpha
    i.particle_radius : FloatSocket
        Particle Radius
    i.particle_alpha : FloatSocket
        Particle alpha
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
    o.normal : VectorSocket
        Normal of the face the current point is nearest
    """

    _name = "Simulate on Faces"
    _asset_name = "Simulate on Faces"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        points: GeometrySocket
        """Points"""
        selection: BooleanSocket
        """Selection"""
        substeps: IntegerSocket
        """Substeps"""
        force: VectorSocket
        """Force"""
        drag: FloatSocket
        """Drag"""
        faces: GeometrySocket
        """Faces"""
        alpha: FloatSocket
        particle_radius: FloatSocket
        """Particle Radius"""
        particle_alpha: FloatSocket
        """Particle alpha"""
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
        normal: VectorSocket
        """Normal of the face the current point is nearest"""

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        points: InputGeometry = None,
        selection: InputBoolean = True,
        substeps: InputInteger = 5,
        force: InputVector = None,
        drag: InputFloat = 0.1,
        faces: InputGeometry = None,
        alpha: InputFloat = 0.0,
        particle_radius: InputFloat = 0.0,
        particle_alpha: InputFloat = 0.0,
        hook_selection: InputBoolean = False,
        hook_target: InputVector = None,
        hook_decay: InputFloat = 2.0,
        pin_selection: InputBoolean = False,
        pin_target: InputVector = None,
    ):
        super().__init__(
            **{
                "Points": points,
                "Selection": selection,
                "Substeps": substeps,
                "Force": force,
                "Drag": drag,
                "Faces": faces,
                "alpha": alpha,
                "Particle Radius": particle_radius,
                "Particle alpha": particle_alpha,
                "Hook Selection": hook_selection,
                "Hook Target": hook_target,
                "Hook Decay": hook_decay,
                "Pin Selection": pin_selection,
                "Pin Target": pin_target,
            }
        )

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        points = tree.inputs.geometry("Points")
        selection = tree.inputs.boolean("Selection", True, hide_value=True)
        substeps = tree.inputs.integer("Substeps", 5, min_value=1, max_value=100)
        with tree.inputs.panel("Simulation"):
            force = tree.inputs.vector(
                "Force",
                (0.0, 0.0, 0.0),
                min_value=-10_000.0,
                max_value=10_000.0,
                hide_value=True,
            )
            drag = tree.inputs.float(
                "Drag", 0.1, min_value=-10_000.0, max_value=10_000.0
            )
            with tree.inputs.panel("Faces"):
                faces = tree.inputs.geometry("Faces")
                alpha = tree.inputs.float(
                    "alpha", 0.0, min_value=-10_000.0, max_value=10_000.0
                )
            with tree.inputs.panel("Particle"):
                particle_radius = tree.inputs.float(
                    "Particle Radius", 0.0, min_value=0.0
                )
                particle_alpha = tree.inputs.float("Particle alpha", 0.0, min_value=0.0)
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
                "Hook Decay", 2.0, min_value=0.0, max_value=10_000.0
            )
        with tree.inputs.panel("Pin"):
            pin_selection = tree.inputs.boolean("Pin Selection", False, hide_value=True)
            pin_target = tree.inputs.vector(
                "Pin Target", (0.0, 0.0, 0.0), hide_value=True, default_input="POSITION"
            )
        geometry = tree.outputs.geometry("Geometry")
        normal = tree.outputs.vector(
            "Normal", description="Normal of the face the current point is nearest"
        )

        boolean_math = g.BooleanMath.subtract(selection, pin_selection)
        store_named_attribute = g.StoreNamedAttribute.point.float(
            points, ~g.NamedAttribute.float("mass").o.exists, "mass", 1.0
        )
        simulation_zone = g.SimulationZone()
        geometry_1 = simulation_zone.items.geometry("Geometry", store_named_attribute)
        math_1 = simulation_zone.delta_time / substeps
        store_named_attribute_1 = g.StoreNamedAttribute.point.float(
            geometry_1.current, name="inverse_mass", value=1.0 / Mass()
        )
        repeat_zone = g.RepeatZone(substeps)
        geometry_2 = repeat_zone.items.geometry("Geometry", store_named_attribute_1)
        group = XPBDInit(
            geometry=geometry_2.current,
            selection=boolean_math,
            force=force,
            drag=drag,
            deltat=math_1,
        )
        group_1 = XPBDSolvePoints(
            geometry=XPBDSolveOnFaces(
                Geometry=group, Faces=faces, alpha=alpha, deltaT=math_1
            ),
            selection=boolean_math,
            radius=particle_radius,
            alpha=particle_alpha,
            deltat=math_1,
        )
        set_position = XPBDSolveHook(
            geometry=group_1,
            selection=hook_selection,
            target=hook_target,
            decay=hook_decay,
            deltat=math_1,
        ) >> g.SetPosition(selection=pin_selection, position=pin_target)
        (
            XPBDFinalise(geometry=set_position, selection=boolean_math, deltat=math_1)
            >> geometry_2.next
        )
        geometry_2.result >> geometry_1.next
        capture = g.CaptureAttribute.point(geometry=geometry_1.result)
        normal_1 = capture.items.vector(
            "Normal", g.SampleNearestSurface.vector(faces, g.Normal().o.normal).o.value
        )

        capture.o.geometry >> geometry
        normal_1.output >> normal


ASSET = SimulateOnFaces

ASSET_METADATA = {
    "catalog_id": "c2c958af-5095-4fc2-884d-709bba965fc4",
}
