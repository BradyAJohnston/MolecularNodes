"""
Importing molecular dynamics trajectories and associated files.
"""

__name__ = "MolecularNodes.trajectory"
__author__ = "Brady Johnston"

import bpy
import MDAnalysis as mda

from ... import blender as bl
from ...session import MNSession
from .trajectory import Trajectory
from . import dna
from bpy.props import StringProperty

bpy.types.Scene.MN_import_md_topology = StringProperty(
    name="Topology",
    description="File path for the toplogy file for the trajectory",
    subtype="FILE_PATH",
    maxlen=0,
)
bpy.types.Scene.MN_import_md_trajectory = StringProperty(
    name="Trajectory",
    description="File path for the trajectory file for the trajectory",
    subtype="FILE_PATH",
    maxlen=0,
)
bpy.types.Scene.MN_import_md_name = StringProperty(
    name="Name",
    description="Name of the molecule on import",
    default="NewTrajectory",
    maxlen=0,
)


def load(top, traj, name="NewTrajectory", style="spheres"):
    top = bl.path_resolve(top)
    traj = bl.path_resolve(traj)

    universe = mda.Universe(top, traj)

    traj = Trajectory(universe=universe)

    traj.create_object(name=name, style=style)

    return traj


class MN_OT_Reload_Trajectory(bpy.types.Operator):
    bl_idname = "mn.reload_trajectory"
    bl_label = "Reload Trajectory"
    bl_description = (
        "Reload the `mda.UNiverse` of the current Object to renable updating"
    )
    bl_options = {"REGISTER", "UNDO"}

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        traj = context.scene.MNSession.match(obj)
        return not traj

    def execute(self, context):
        obj = context.active_object
        session: MNSession = context.scene.MNSession
        topo = obj.mn.filepath_topology
        traj = obj.mn.filepath_trajectory

        if "oxdna" in obj.mn.entity_type:
            uni = mda.Universe(
                topo, traj, topology_format=dna.OXDNAParser, format=dna.OXDNAReader
            )
            traj = dna.OXDNA(uni)
        else:
            uni = mda.Universe(topo, traj)
            traj = Trajectory(uni)

        traj.object = obj
        traj.set_frame(context.scene.frame_current)
        return {"FINISHED"}


class MN_OT_Import_Trajectory(bpy.types.Operator):
    bl_idname = "mn.import_trajectory"
    bl_label = "Import Protein MD"
    bl_description = "Load molecular dynamics trajectory"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        scene = context.scene
        topology_file = scene.MN_import_md_topology
        trajectory_file = scene.MN_import_md_trajectory
        name = scene.MN_import_md_name

        trajectory = load(
            top=topology_file,
            traj=trajectory_file,
            name=name,
            style=scene.mn.import_style,
        )

        context.view_layer.objects.active = trajectory.object
        context.scene.frame_start = 0
        context.scene.frame_end = trajectory.universe.trajectory.n_frames

        self.report(
            {"INFO"},
            message=f"Imported '{topology_file}' as {trajectory.object.name} "
            f"with {str(trajectory.universe.trajectory.n_frames)} "
            f"frames from '{trajectory_file}'.",
        )

        return {"FINISHED"}


# UI


def panel(layout, scene):
    layout.label(text="Load MD Trajectories", icon="FILE_TICK")
    layout.separator()
    col = layout.column(align=True)
    row_import = col.row()
    row_import.prop(scene, "MN_import_md_name")
    row_import.operator("mn.import_trajectory", text="Load")
    col.separator()
    col.prop(scene, "MN_import_md_topology")
    col.prop(scene, "MN_import_md_trajectory")

    layout.separator()
    layout.label(text="Options", icon="MODIFIER")
    row = layout.row()
    row.prop(scene.mn, "import_node_setup", text="")
    col = row.column()
    col.prop(scene.mn, "import_style")
    col.enabled = scene.mn.import_node_setup


CLASSES = [MN_OT_Import_Trajectory, MN_OT_Reload_Trajectory]
