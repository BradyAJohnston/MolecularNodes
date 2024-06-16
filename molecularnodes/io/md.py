"""
Importing molecular dynamics trajectories and associated files.
"""

__name__ = "MolecularNodes.trajectory"
__author__ = "Brady Johnston"

import bpy
import MDAnalysis as mda

from ..blender import path_resolve
from .parse.mda import MNUniverse
from bpy.props import StringProperty, IntProperty, BoolProperty

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
bpy.types.Scene.MN_import_md_frame_start = IntProperty(
    name="Start", description="Frame start for importing MD trajectory", default=0
)
bpy.types.Scene.MN_import_md_frame_step = IntProperty(
    name="Step", description="Frame step for importing MD trajectory", default=1
)
bpy.types.Scene.MN_import_md_frame_stop = IntProperty(
    name="Stop", description="Frame stop for importing MD trajectory", default=499
)
bpy.types.Scene.MN_md_selection = StringProperty(
    name="Import Filter",
    description='Custom MDAnalysis selection string, removing unselecte atoms. See: "https://docs.mdanalysis.org/stable/documentation_pages/selections.html"',
    default="all",
)
bpy.types.Scene.MN_md_in_memory = BoolProperty(
    name="In Memory",
    description="True will load all of the requested frames into the scene and into memory. False will stream the trajectory from a live MDAnalysis session",
    default=False,
)
bpy.types.Scene.list_index = IntProperty(
    name="Index for trajectory selection list.", default=0
)


def load(
    top,
    traj,
    name="NewTrajectory",
    style="spheres",
    start: int = 0,
    step: int = 1,
    stop: int = 499,
    subframes: int = 0,
    in_memory: bool = False,
):
    top = path_resolve(top)
    traj = path_resolve(traj)

    universe = mda.Universe(top, traj)

    mn_universe = MNUniverse(universe=universe)

    mn_universe.create_model(name=name, style=style, subframes=subframes)
    bpy.context.scene.MNSession.universes.append(mn_universe)

    return mn_universe


class MN_OT_Import_Protein_MD(bpy.types.Operator):
    bl_idname = "mn.import_protein_md"
    bl_label = "Import Protein MD"
    bl_description = "Load molecular dynamics trajectory"
    bl_options = {"REGISTER", "UNDO"}

    @classmethod
    def poll(cls, context):
        return True

    def execute(self, context):
        scene = context.scene
        top = scene.MN_import_md_topology
        traj = scene.MN_import_md_trajectory
        name = scene.MN_import_md_name

        mu = load(
            top=top,
            traj=traj,
            name=name,
            style=scene.MN_import_style,
            start=scene.MN_import_md_frame_start,
            stop=scene.MN_import_md_frame_stop,
            step=scene.MN_import_md_frame_step,
            in_memory=scene.MN_md_in_memory,
        )

        bpy.context.view_layer.objects.active = mu.object

        self.report(
            {"INFO"},
            message=f"Imported '{top}' as {name} "
            f"with {str(mu.universe.trajectory.n_frames)} "
            f"frames from '{traj}'.",
        )

        return {"FINISHED"}


# UI


class TrajectorySelectionItem(bpy.types.PropertyGroup):
    """Group of properties for custom selections for MDAnalysis import."""

    name: StringProperty(  # type: ignore
        name="Name",
        description="Becomes the attribute name when applied to the mesh",
        default="custom_selection",
    )

    text: StringProperty(  # type: ignore
        name="Selection String",
        description="String that provides a selection through MDAnalysis",
        default="name CA",
    )

    update: BoolProperty(  # type: ignore
        name="Update",
        description="Recalculate the selection on frame change",
        default=True,
    )

    valid: BoolProperty(  # type: ignore
        name="Valid",
        description="If the previous attempt to calculate the selection succeeded",
        default=True,
    )

    periodic: BoolProperty(  # type: ignore
        name="Periodic",
        description="For geometric selections, whether to account for atoms in different periodic images when searching",
        default=True,
    )


class MN_UL_TrajectorySelectionListUI(bpy.types.UIList):
    """UI List"""

    def draw_item(
        self, context, layout, data, item, icon, active_data, active_propname, index
    ):
        custom_icon = "VIS_SEL_11"

        if self.layout_type in {"DEFAULT", "COMPACT"}:
            if not item.valid:
                custom_icon = "ERROR"
            layout.label(text=item.name, icon=custom_icon)

        elif self.layout_type in {"GRID"}:
            layout.alignment = "CENTER"
            layout.label(text="", icon=custom_icon)


class TrajectorySelection_OT_NewItem(bpy.types.Operator):
    """Add a new custom selection to the list."""

    bl_idname = "trajectory_selection_list.new_item"
    bl_label = "+"

    def execute(self, context):
        context.scene.trajectory_selection_list.add()
        return {"FINISHED"}


class MDASelection_OT_NewItem(bpy.types.Operator):
    "Add a new custom selection to a trajectory"

    bl_idname = "mda.new_item"
    bl_label = "+"

    def execute(self, context):
        o = context.active_object
        o.mda.add()
        o["list_index"] = len(o.mda) - 1

        return {"FINISHED"}


class MDASelection_OT_DeleteItem(bpy.types.Operator):
    bl_idname = "mda.delete_item"
    bl_label = "-"

    @classmethod
    def poll(cls, context):
        return context.active_object.mda

    def execute(self, context):
        o = context.active_object
        my_list = o.mda
        index = context.scene.list_index

        my_list.remove(index)
        context.scene.list_index = len(my_list) - 1

        return {"FINISHED"}


class TrajectorySelection_OT_DeleteIem(bpy.types.Operator):
    bl_idname = "trajectory_selection_list.delete_item"
    bl_label = "-"

    @classmethod
    def poll(cls, context):
        return context.scene.trajectory_selection_list

    def execute(self, context):
        my_list = context.scene.trajectory_selection_list
        index = context.scene.list_index

        my_list.remove(index)
        context.scene.list_index = min(max(0, index - 1), len(my_list) - 1)

        return {"FINISHED"}


# def custom_selections(layout, scene):
#     layout.label(text="Custom Selections")
#     row = layout.row(align=True)

#     row = row.split(factor=0.9)
#     row.template_list(
#         "MN_UL_TrajectorySelectionListUI",
#         "A list",
#         scene,
#         "trajectory_selection_list",
#         scene,
#         "list_index",
#         rows=3,
#     )
# col = row.column()
# col.operator("trajectory_selection_list.new_item", icon="ADD", text="")
# col.operator("trajectory_selection_list.delete_item", icon="REMOVE", text="")
# if scene.list_index >= 0 and scene.trajectory_selection_list:
#     item = scene.trajectory_selection_list[scene.list_index]

#     col = layout.column(align=False)
#     col.separator()

# col.prop(item, "name")
# col.prop(item, "selection")
# col.prop(item, "update")


def panel(layout, scene):
    layout.label(text="Load MD Trajectories", icon="FILE_TICK")
    layout.separator()
    col = layout.column(align=True)
    row_import = col.row()
    row_import.prop(scene, "MN_import_md_name")
    row_import.operator("mn.import_protein_md", text="Load")
    col.separator()
    col.prop(scene, "MN_import_md_topology")
    col.prop(scene, "MN_import_md_trajectory")

    layout.separator()
    layout.label(text="Options", icon="MODIFIER")
    row = layout.row()
    row.prop(scene, "MN_import_node_setup", text="")
    col = row.column()
    col.prop(scene, "MN_import_style")
    col.enabled = scene.MN_import_node_setup
    layout.prop(scene, "MN_md_selection")
    row_frame = layout.row(heading="Frames", align=True)
    row_frame.prop(scene, "MN_md_in_memory")
    row = row_frame.row(align=True)
    row.prop(scene, "MN_import_md_frame_start")
    row.prop(scene, "MN_import_md_frame_step")
    row.prop(scene, "MN_import_md_frame_stop")
    row.enabled = scene.MN_md_in_memory
    # custom_selections(layout, scene)


CLASSES = [
    # TrajectorySelectionItem,  # has to be registered before the others to work properly
    # MN_UL_TrajectorySelectionListUI,
    # TrajectorySelection_OT_DeleteIem,
    # TrajectorySelection_OT_NewItem,
    MN_OT_Import_Protein_MD,
]
