from MDAnalysis import Universe
import bpy
from ... import color
from ...blender import coll, nodes
from ... import bpyd
from ...bpyd import AttributeTypes

from .oxdna.OXDNAParser import OXDNAParser
from .oxdna.OXDNAReader import OXDNAReader
from .trajectory import Trajectory

DNA_SCALE = 10

bpy.types.Scene.MN_import_oxdna_topology = bpy.props.StringProperty(
    name="Toplogy",
    description="File path for the topology to import (.top)",
    subtype="FILE_PATH",
    maxlen=0,
)
bpy.types.Scene.MN_import_oxdna_trajectory = bpy.props.StringProperty(
    name="Trajectory",
    description="File path for the trajectory to import (.oxdna / .dat)",
    subtype="FILE_PATH",
    maxlen=0,
)
bpy.types.Scene.MN_import_oxdna_name = bpy.props.StringProperty(
    name="Name",
    description="Name of the created object.",
    default="NewOrigami",
    maxlen=0,
)


class OXDNA(Trajectory):
    def __init__(self, universe: Universe, world_scale: float = 0.01):
        super().__init__(universe=universe, world_scale=world_scale)
        self._att_names = (
            "base_vector",
            "base_normal",
            "velocity",
            "angular_velocity",
        )

    def create_object(
        self,
        style: str = "oxdna",
        name: str = "NewUniverseObject",
    ):
        self.object = bpyd.create_object(
            name=name,
            collection=coll.mn(),
            vertices=self.univ_positions * self.world_scale * DNA_SCALE,
            edges=self.bonds,
        )
        self.object.mn.uuid = self.uuid
        self.object.mn.molecule_type = "md"
        self._update_timestep_values()

        for name in ("chain_id", "res_num", "res_id"):
            self.store_named_attribute(getattr(self, name), name)

        self.store_named_attribute(
            data=color.color_chains_equidistant(self.chain_id),
            name="Color",
            atype=AttributeTypes.FLOAT_COLOR,
        )

        if style:
            nodes.create_starting_node_tree(self.object, style="oxdna", color=None)

        return self.object

    def _update_positions(self, frame: int) -> None:
        super()._update_positions(frame)
        self._update_timestep_values()

    def _update_timestep_values(self):
        for name in self._att_names:
            try:
                self.store_named_attribute(
                    self.universe.trajectory.ts.data[name] * self.world_scale, name=name
                )
            except KeyError:
                pass


def load(top, traj, name="oxDNA", style="oxdna", world_scale=0.01):
    univ = Universe(top, traj, topology_format=OXDNAParser, format=OXDNAReader)
    traj = OXDNA(univ, world_scale=world_scale * DNA_SCALE)
    traj.create_object(name=name, style=style)
    return traj


class MN_OT_Import_OxDNA_Trajectory(bpy.types.Operator):
    bl_idname = "mn.import_oxdna"
    bl_label = "Load"
    bl_description = "Will import the given file and toplogy."
    bl_options = {"REGISTER"}

    def execute(self, context):
        s = context.scene
        load(
            top=s.MN_import_oxdna_topology,
            traj=s.MN_import_oxdna_trajectory,
            name=s.MN_import_oxdna_name,
        )
        return {"FINISHED"}


def panel(layout, scene):
    layout.label(text="Load oxDNA File", icon="FILE_TICK")
    layout.separator()
    row = layout.row()
    row.prop(scene, "MN_import_oxdna_name")
    row.operator("mn.import_oxdna")
    col = layout.column(align=True)
    col.prop(scene, "MN_import_oxdna_topology")
    col.prop(scene, "MN_import_oxdna_trajectory")
