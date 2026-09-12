# Node-group asset "Simulate Elastic Network" (GeometryNodeTree), dumped by nodebpy.assets.dump_library.
# Rebuild the library with nodebpy.assets.build_library (python -m nodebpy.assets build).
# _build_group() is the source of truth: the docstring, __init__ and accessors are regenerated from it on the next dump.
from typing import TYPE_CHECKING, Literal
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
    MenuSocket,
    PackageLibrary,
    SocketAccessor,
    VectorSocket,
)
from nodebpy.types import (
    InputBoolean,
    InputFloat,
    InputGeometry,
    InputInteger,
    InputMenu,
    InputVector,
)
from ._shared.constraint_distance import ConstraintDistance
from ._shared.inverse_mass import InverseMass
from ._shared.xpbd_finalise import XPBDFinalise
from ._shared.xpbd_init import XPBDInit
from ._shared.xpbd_solve_hook import XPBDSolveHook
from ._shared.xpbd_solve_points import XPBDSolvePoints
from .edge_info import EdgeInfo
from .edge_length import EdgeLength
from .mass import Mass


class XPBDSolveEdges(CustomGeometryGroup):
    _name = "XPBD Solve Edges"
    _color_tag = "GEOMETRY"

    def _build_group(self, tree: TreeBuilder[GeometryNodeTree]) -> None:
        geometry = tree.inputs.geometry("Geometry")
        selection = tree.inputs.boolean("Selection", True, hide_value=True)
        distance = tree.inputs.float(
            "Distance", 0.27, min_value=0.0, max_value=10_000.0
        )
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

        group = InverseMass()
        repeat_zone = g.RepeatZone(
            g.AttributeStatistic.point.float(
                geometry, attribute=g.EdgesOfVertex().o.total
            ).o.max
        )
        geometry_2 = repeat_zone.items.geometry("Geometry", geometry)
        correction = repeat_zone.items.vector("Correction")
        value = repeat_zone.items.integer("Value")
        group_1 = EdgeInfo(edge_index=repeat_zone.iteration)
        boolean_math = selection & group_1.o.is_valid
        group_2 = ConstraintDistance(
            target=group_1.o.point_position,
            distance=distance.edge.at(group_1.o.edge_index),
            w1=group.o.w,
            w2=group.o.w.point.at(group_1.o.edge_index),
            alpha=alpha.edge.at(group_1.o.edge_index),
            deltat=deltat,
        )
        geometry_2.current >> geometry_2.next
        (
            correction.current
            + boolean_math.switch.vector((0.0, 0.0, 0.0), group_2.o.correction)
            >> correction.next
        )
        g.IntegerMath(value=value.current, value_001=boolean_math) >> value.next
        (
            geometry_2.result
            >> g.SetPosition(offset=correction.result / value.result)
            >> geometry_1
        )


class SimulateElasticNetwork(AssetGeometryGroup):
    """
    Simulate Elastic Network

    Parameters
    ----------
    mesh : InputGeometry
        Mesh
    selection : InputBoolean
        Selection
    substeps : InputInteger
        Substeps
    force : InputVector
        Force
    drag : InputFloat
        Drag
    edge_length_source : InputMenu | Literal["Original", "Custom"]
        Edge Length Source
    edge_length : InputFloat
        Edge Length
    edge_alpha : InputFloat
        Edge alpha
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
    i.mesh : GeometrySocket
        Mesh
    i.selection : BooleanSocket
        Selection
    i.substeps : IntegerSocket
        Substeps
    i.force : VectorSocket
        Force
    i.drag : FloatSocket
        Drag
    i.edge_length_source : MenuSocket
        Edge Length Source
    i.edge_length : FloatSocket
        Edge Length
    i.edge_alpha : FloatSocket
        Edge alpha
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
    """

    _name = "Simulate Elastic Network"
    _asset_name = "Simulate Elastic Network"
    _library = PackageLibrary(__file__, "../../assets/nodes.blend")
    _color_tag = "GEOMETRY"

    class _Inputs(SocketAccessor):
        mesh: GeometrySocket
        """Mesh"""
        selection: BooleanSocket
        """Selection"""
        substeps: IntegerSocket
        """Substeps"""
        force: VectorSocket
        """Force"""
        drag: FloatSocket
        """Drag"""
        edge_length_source: MenuSocket
        """Edge Length Source"""
        edge_length: FloatSocket
        """Edge Length"""
        edge_alpha: FloatSocket
        """Edge alpha"""
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

    if TYPE_CHECKING:

        @property
        def i(self) -> _Inputs: ...
        @property
        def o(self) -> _Outputs: ...

    def __init__(
        self,
        mesh: InputGeometry = None,
        selection: InputBoolean = True,
        substeps: InputInteger = 5,
        force: InputVector = None,
        drag: InputFloat = 0.1,
        edge_length_source: InputMenu | Literal["Original", "Custom"] = "Original",
        edge_length: InputFloat = 0.1,
        edge_alpha: InputFloat = 0.0,
        particle_radius: InputFloat = 0.005,
        particle_alpha: InputFloat = 0.0,
        hook_selection: InputBoolean = False,
        hook_target: InputVector = None,
        hook_decay: InputFloat = 2.0,
        pin_selection: InputBoolean = False,
        pin_target: InputVector = None,
    ):
        super().__init__(
            **{
                "Mesh": mesh,
                "Selection": selection,
                "Substeps": substeps,
                "Force": force,
                "Drag": drag,
                "Edge Length Source": edge_length_source,
                "Edge Length": edge_length,
                "Edge alpha": edge_alpha,
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
        mesh = tree.inputs.geometry("Mesh")
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
            drag = tree.inputs.float("Drag", 0.1, min_value=0.0, max_value=10_000.0)
            with tree.inputs.panel("Edge"):
                edge_length_source = tree.inputs.menu(
                    "Edge Length Source", optional_label=True
                )
                edge_length = tree.inputs.float(
                    "Edge Length", 0.1, min_value=0.0, max_value=10_000.0
                )
                edge_alpha = tree.inputs.float("Edge alpha", 0.0, min_value=0.0)
            with tree.inputs.panel("Particle"):
                particle_radius = tree.inputs.float(
                    "Particle Radius", 0.005, min_value=0.0
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

        boolean_math = g.BooleanMath.subtract(selection, pin_selection)
        menu_switch = g.MenuSwitch.float(
            edge_length_source,
            {
                "Original": g.NamedAttribute.float("tmp_length").o.attribute,
                "Custom": edge_length,
            },
        )
        store_named_attribute = (
            mesh
            >> g.StoreNamedAttribute.point.float(
                selection=~g.NamedAttribute.float("mass").o.exists,
                name="mass",
                value=1.0,
            )
            >> g.StoreNamedAttribute.edge.float(name="tmp_length", value=EdgeLength())
        )
        simulation_zone = g.SimulationZone()
        geometry_1 = simulation_zone.items.geometry("Geometry", store_named_attribute)
        math_1 = simulation_zone.delta_time / substeps
        store_named_attribute_1 = g.StoreNamedAttribute.point.float(
            geometry_1.current, name="inverse_mass", value=Mass().o.mass / 0.5
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
        group_1 = XPBDSolveEdges(
            Geometry=group,
            Selection=boolean_math,
            Distance=menu_switch.o.output,
            alpha=edge_alpha,
            deltaT=math_1,
        )
        group_2 = XPBDSolvePoints(
            geometry=group_1,
            selection=boolean_math,
            radius=particle_radius,
            alpha=particle_alpha,
            deltat=math_1,
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
            >> geometry_2.next
        )
        geometry_2.result >> geometry_1.next
        (
            geometry_1.result
            >> g.RemoveNamedAttribute(pattern_mode="Wildcard", name="tmp_*")
            >> geometry
        )

        edge_length_source.default_value = "Original"


ASSET = SimulateElasticNetwork

ASSET_METADATA = {
    "catalog_id": "c2c958af-5095-4fc2-884d-709bba965fc4",
}
