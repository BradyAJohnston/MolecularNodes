"""
Importing molecular dynamics trajectories and associated files.
"""

__name__ = "MolecularNodes.trajectory"
__author__ = "Brady Johnston"

import bpy
import MDAnalysis as mda

from ..blender import path_resolve
from .parse.mda import MNUniverse, _update_universes
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


def load(
    top,
    traj,
    name="NewTrajectory",
    style="spheres",
    subframes: int = 0,
):
    top = path_resolve(top)
    traj = path_resolve(traj)

    universe = mda.Universe(top, traj)

    mn_universe = MNUniverse(universe=universe)

    mn_universe.create_model(name=name, style=style, subframes=subframes)

    return mn_universe


class MN_OT_Import_Protein_MD(bpy.types.Operator):
    bl_idname = "mn.import_protein_md"
    bl_label = "Import Protein MD"
    bl_description = "Load molecular dynamics trajectory"
    bl_options = {"REGISTER", "UNDO"}

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
        )

        context.view_layer.objects.active = mu.object
        context.scene.frame_start = 0
        context.scene.frame_end = mu.universe.trajectory.n_frames

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
        description="Name of the attribute on the mesh",
        default="custom_selection",
        update=_update_universes,
    )

    selection_str: StringProperty(  # type: ignore
        name="Selection",
        description="Selection to be applied, written in the MDAnalysis selection language",
        default="name CA",
        update=_update_universes,
    )

    updating: BoolProperty(  # type: ignore
        name="Updating",
        description="Recalculate the selection on scene frame change",
        default=True,
        update=_update_universes,
    )

    periodic: BoolProperty(  # type: ignore
        name="Periodic",
        description="For geometric selections, whether to account for atoms in different periodic images when searching",
        default=True,
        update=_update_universes,
    )
    message: StringProperty(  # type: ignore
        name="Message",
        description="Message to report back from `universe.select_atoms()`",
        default="",
    )


class MN_UL_TrajectorySelectionListUI(bpy.types.UIList):
    """UI List"""

    def draw_item(
        self, context, layout, data, item, icon, active_data, active_propname, index
    ):
        custom_icon = "VIS_SEL_11"

        if self.layout_type in {"DEFAULT", "COMPACT"}:
            row = layout.row()
            if item.message != "":
                custom_icon = "ERROR"
                row.alert = True

            row.prop(item, "name", text="", emboss=False)
            row.prop(item, "updating", icon_only=True, icon="FILE_REFRESH")
            row.prop(item, "periodic", icon_only=True, icon="CUBE")

        elif self.layout_type in {"GRID"}:
            layout.alignment = "CENTER"
            layout.label(text="", icon=custom_icon)


class MN_OT_Universe_Selection_Add(bpy.types.Operator):
    "Add a new custom selection to a trajectory"

    bl_idname = "mn.universe_selection_add"
    bl_label = "+"
    bl_description = "Add a new boolean attribute for the given MDA selection string"

    def execute(self, context):
        bob = context.active_object
        bob.mn_universe_selections.add()
        i = int(len(bob.mn_universe_selections) - 1)
        bob.mn_universe_selections[i].name = f"selection_{i + 1}"
        bob.mn["list_index"] = i
        _update_universes(self, context)

        return {"FINISHED"}


class MN_OT_Universe_Selection_Delete(bpy.types.Operator):
    bl_idname = "mda.delete_item"
    bl_label = "-"
    bl_description = "Delete the given boolean selection from the universe"

    @classmethod
    def poll(cls, context):
        return context.active_object.mn_universe_selections

    def execute(self, context):
        bob = context.active_object
        index = bob.mn.universe_selection_index

        sel_list = bob.mn_universe_selections
        sel_list.remove(index)
        bob.mn.universe_selection_index = len(sel_list) - 1
        _update_universes(self, context)

        return {"FINISHED"}


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


CLASSES = [
    TrajectorySelectionItem,  # has to be registered before the others to work properly
    MN_UL_TrajectorySelectionListUI,
    MN_OT_Universe_Selection_Add,
    MN_OT_Universe_Selection_Delete,
    MN_OT_Import_Protein_MD,
]
