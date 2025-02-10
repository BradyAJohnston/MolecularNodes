from MDAnalysis import Universe
import bpy
from bpy.props import StringProperty
from .ops import TrajectoryImportOperator
from ... import color
from ...blender import coll, nodes
import databpy
from databpy import AttributeTypes

from enum import Enum

from ..entity import EntityType

from .oxdna.OXDNAParser import OXDNAParser
from .oxdna.OXDNAReader import OXDNAReader
from .trajectory import Trajectory

DNA_SCALE = 10


class OXDNA(Trajectory):
    def __init__(self, universe: Universe, world_scale: float = 0.01):
        super().__init__(universe=universe, world_scale=world_scale * DNA_SCALE)
        self._entity_type = EntityType.MD_OXDNA
        self._att_names = (
            "base_vector",
            "base_normal",
            "velocity",
            "angular_velocity",
        )

    def _create_object(self, style: str = "oxdna", name: str = "NewUniverseObject"):
        self.object = databpy.create_object(
            name=name,
            collection=coll.mn(),
            vertices=self.univ_positions,
            edges=self.bonds,
        )
        self._update_timestep_values()

        for name in ("chain_id", "res_id", "res_name"):
            if name == "res_name":
                att_name = "res_num"
            else:
                att_name = name
            self.store_named_attribute(getattr(self, att_name), name)

        self.store_named_attribute(
            data=color.color_chains_equidistant(self.chain_id),
            name="Color",
            atype=AttributeTypes.FLOAT_COLOR,
        )

        if style:
            nodes.create_starting_node_tree(self.object, style="oxdna", color=None)

        return self.object

    def set_frame(self, frame: int) -> None:
        super()._update_positions(frame)
        self._update_timestep_values()

    def _update_timestep_values(self):
        for name in self._att_names:
            try:
                self.store_named_attribute(
                    self.universe.trajectory.ts.data[name] * self.world_scale, name=name
                )
            except KeyError as e:
                print(e)


def load(top, traj, name="oxDNA", style="oxdna", world_scale=0.01):
    univ = Universe(top, traj, topology_format=OXDNAParser, format=OXDNAReader)
    traj = OXDNA(univ, world_scale=world_scale)
    traj.create_object(name=name, style=style)
    return traj


class MN_OT_Import_OxDNA_Trajectory(TrajectoryImportOperator):
    bl_idname = "mn.import_oxdna"

    def execute(self, context):
        load(top=self.topology, traj=self.trajectory, name=self.name)
        return {"FINISHED"}


def panel(layout, scene):
    layout.label(text="Load oxDNA File", icon="FILE_TICK")
    layout.separator()
    row = layout.row()
    row.prop(scene.mn, "import_oxdna_name")
    op = row.operator("mn.import_oxdna")
    op.name = scene.mn.import_oxdna_name
    op.topology = scene.mn.import_oxdna_topology
    op.trajectory = scene.mn.import_oxdna_trajectory
    col = layout.column(align=True)
    col.prop(scene.mn, "import_oxdna_topology")
    col.prop(scene.mn, "import_oxdna_trajectory")
