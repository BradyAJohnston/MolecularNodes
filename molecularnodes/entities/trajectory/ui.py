"""
Importing molecular dynamics trajectories and associated files.
"""

__name__ = "MolecularNodes.trajectory"
__author__ = "Brady Johnston"

import bpy
import MDAnalysis as mda

from pathlib import Path
from ... import blender as bl
from ...style import STYLE_ITEMS
from ...session import MNSession
from .base import Trajectory
from . import dna
from bpy.props import StringProperty, EnumProperty, BoolProperty


def load(
    top: str | Path,
    traj: str | Path,
    name: str = "NewTrajectory",
    style: str | None = "spheres",
):
    top = bl.path_resolve(top)
    traj = bl.path_resolve(traj)

    universe = mda.Universe(top, traj)
    trajectory = Trajectory(universe=universe)
    trajectory.create_object(name=name, style=style)

    return trajectory


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

    topology: StringProperty(  # type: ignore
        name="Topology",
        description="File path for the toplogy file for the trajectory",
        subtype="FILE_PATH",
        maxlen=0,
    )
    trajectory: StringProperty(  # type: ignore
        name="Trajectory",
        description="File path for the trajectory file for the trajectory",
        subtype="FILE_PATH",
        maxlen=0,
    )
    name: StringProperty(  # type: ignore
        name="Name",
        description="Name of the molecule on import",
        default="NewTrajectory",
        maxlen=0,
    )
    style: EnumProperty(  # type: ignore
        name="Style",
        description="Default style for importing",
        items=STYLE_ITEMS,
        default="spheres",
    )
    setup_nodes: BoolProperty(  # type: ignore
        name="Setup Nodes",
        description="Add nodes to the scene to load the trajectory",
        default=True,
    )

    def execute(self, context):
        trajectory = load(
            top=self.topology,
            traj=self.trajectory,
            name=self.name,
            style=self.style if self.setup_nodes else None,
        )

        context.view_layer.objects.active = trajectory.object
        context.scene.frame_start = 0
        context.scene.frame_end = trajectory.universe.trajectory.n_frames

        self.report(
            {"INFO"},
            message=f"Imported '{self.topology}' as {trajectory.name} "
            f"with {str(trajectory.universe.trajectory.n_frames)} "
            f"frames from '{self.trajectory}'.",
        )

        return {"FINISHED"}


# UI


def panel(layout, scene):
    layout.label(text="Load MD Trajectories", icon="FILE_TICK")
    layout.separator()
    col = layout.column(align=True)
    row_import = col.row()
    row_import.prop(scene.mn, "import_md_name")
    op = row_import.operator("mn.import_trajectory", text="Load")
    op.topology = scene.mn.import_md_topology
    op.trajectory = scene.mn.import_md_trajectory
    op.name = scene.mn.import_md_name
    op.style = scene.mn.import_style
    op.setup_nodes = scene.mn.import_node_setup
    col.separator()
    col.prop(scene.mn, "import_md_topology")
    col.prop(scene.mn, "import_md_trajectory")

    layout.separator()
    layout.label(text="Options", icon="MODIFIER")
    row = layout.row()
    row.prop(scene.mn, "import_node_setup", text="")
    col = row.column()
    col.prop(scene.mn, "import_style")
    col.enabled = scene.mn.import_node_setup


CLASSES = [MN_OT_Import_Trajectory, MN_OT_Reload_Trajectory]
