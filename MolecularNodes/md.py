"""
Importing molecular dynamics trajectories and associated files.
"""

__name__ = "MolecularNodes.trajectory"
__author__ = "Brady Johnston"

import bpy
import numpy as np
import warnings
import MDAnalysis as mda

from . import data
from . import coll
from . import obj
from . import nodes
from .mda import MDAnalysisSession

bpy.types.Scene.mol_import_md_topology = bpy.props.StringProperty(
    name = 'path_topology', 
    description = 'File path for the toplogy file for the trajectory', 
    options = {'TEXTEDIT_UPDATE'}, 
    default = '/nethome/yzhuang/jupyter_playground/blender/md.tpr', 
    subtype = 'FILE_PATH', 
    maxlen = 0
    )
bpy.types.Scene.mol_import_md_trajectory = bpy.props.StringProperty(
    name = 'path_trajectory', 
    description = 'File path for the trajectory file for the trajectory', 
    options = {'TEXTEDIT_UPDATE'}, 
    default = '/nethome/yzhuang/jupyter_playground/blender/md.xtc', 
    subtype = 'FILE_PATH', 
    maxlen = 0
    )
bpy.types.Scene.mol_import_md_name = bpy.props.StringProperty(
    name = 'mol_md_name', 
    description = 'Name of the molecule on import', 
    options = {'TEXTEDIT_UPDATE'}, 
    default = 'NewTrajectory', 
    subtype = 'NONE', 
    maxlen = 0
    )
bpy.types.Scene.mol_import_md_frame_start = bpy.props.IntProperty(
    name = "mol_import_md_frame_start", 
    description = "Frame start for importing MD trajectory", 
    subtype = 'NONE',
    default = 0
)
bpy.types.Scene.mol_import_md_frame_step = bpy.props.IntProperty(
    name = "mol_import_md_frame_step", 
    description = "Frame step for importing MD trajectory", 
    subtype = 'NONE',
    default = 1
)
bpy.types.Scene.mol_import_md_frame_end = bpy.props.IntProperty(
    name = "mol_import_md_frame_end", 
    description = "Frame end for importing MD trajectory", 
    subtype = 'NONE',
    default = 49
)
bpy.types.Scene.mol_md_selection = bpy.props.StringProperty(
    name = 'md_selection', 
    description = 'Custom selection string when importing MD simulation. See: "https://docs.mdanalysis.org/stable/documentation_pages/selections.html"', 
    options = {'TEXTEDIT_UPDATE'}, 
    default = 'not (name H* or name OW)', 
    subtype = 'NONE'
    )
bpy.types.Scene.list_index = bpy.props.IntProperty(
    name = "Index for trajectory selection list.", 
    default = 0
)
    

class MOL_OT_Import_Protein_MD(bpy.types.Operator):
    bl_idname = "mol.import_protein_md"
    bl_label = "Import Protein MD"
    bl_description = "Load molecular dynamics trajectory"
    bl_options = {"REGISTER", "UNDO"}

    @classmethod
    def poll(cls, context):
        return True

    def execute(self, context):
        file_top = bpy.context.scene.mol_import_md_topology
        file_traj = bpy.context.scene.mol_import_md_trajectory
        name = bpy.context.scene.mol_import_md_name
        selection = bpy.context.scene.mol_md_selection
        md_start = bpy.context.scene.mol_import_md_frame_start
        md_step =  bpy.context.scene.mol_import_md_frame_step
        md_end =   bpy.context.scene.mol_import_md_frame_end
        include_bonds = bpy.context.scene.mol_import_include_bonds
        custom_selections = bpy.context.scene.trajectory_selection_list
        
        universe = mda.Universe(file_top, file_traj)
        
        self.mda_session = MDAnalysisSession(
            universe = universe, 
            name = name, 
            selection = selection, 
            md_start = md_start, 
            md_step = md_step, 
            md_end = md_end, 
            include_bonds = include_bonds, 
            custom_selections = custom_selections
        )

        bpy.context.view_layer.objects.active = self.mda_session.mol_objects[0]
        self.report(
            {'INFO'}, 
            message=f"Imported '{file_top}' as {self.mda_session.mol_objects[0].name} "
                    f"with {str(self.mda_session.n_frames)} "
                    f"frames from '{file_traj}'."
                )
        
        return {"FINISHED"}


#### UI

class TrajectorySelectionItem(bpy.types.PropertyGroup):
    """Group of properties for custom selections for MDAnalysis import."""
    bl_idname = "testing"
    
    name: bpy.props.StringProperty(
        name="Attribute Name", 
        description="Attribute", 
        default="custom_selection"
    )
    
    selection: bpy.props.StringProperty(
        name="Selection String", 
        description="String that provides a selection through MDAnalysis", 
        default = "name CA"
    )


# have to manually register this class otherwise the PropertyGroup registration fails
bpy.utils.register_class(TrajectorySelectionItem)
bpy.types.Scene.trajectory_selection_list = bpy.props.CollectionProperty(
    type = TrajectorySelectionItem
)

class MOL_UL_TrajectorySelectionListUI(bpy.types.UIList):
    """UI List"""
    
    def draw_item(self, context, layout, data, item, 
                  icon, active_data, active_propname, index):
        custom_icon = "VIS_SEL_11"
        
        if self.layout_type in {'DEFAULT', 'COMPACT'}:
            layout.label(text = item.name, icon = custom_icon)
        
        elif self.layout_type in {'GRID'}:
            layout.alignment = 'CENTER'
            layout.label(text = "", icon = custom_icon)
            

class TrajectorySelection_OT_NewItem(bpy.types.Operator):
    """Add a new custom selection to the list."""
    
    bl_idname = "trajectory_selection_list.new_item"
    bl_label = "+"
    
    def execute(self, context):
        context.scene.trajectory_selection_list.add()
        return {'FINISHED'}

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
        
        return {'FINISHED'}

def panel(layout_function, scene):
    col_main = layout_function.column(heading = '', align = False)
    col_main.alert = False
    col_main.enabled = True
    col_main.active = True
    col_main.label(text = "Import Molecular Dynamics Trajectories")
    row_import = col_main.row()
    row_import.prop(
        bpy.context.scene, 'mol_import_md_name', 
        text = "Name", 
        emboss = True
    )
    row_import.operator('mol.import_protein_md', text = "Load", icon='FILE_TICK')
    row_topology = col_main.row(align = True)
    row_topology.prop(
        bpy.context.scene, 'mol_import_md_topology', 
        text = 'Topology',
        emboss = True
    )
    row_trajectory = col_main.row()
    row_trajectory.prop(
        bpy.context.scene, 'mol_import_md_trajectory', 
        text = 'Trajectory', 
        icon_value = 0, 
        emboss = True
    )
    row_frame = col_main.row(heading = "Frames", align = True)
    row_frame.prop(
        bpy.context.scene, 'mol_import_md_frame_start', 
        text = 'Start',
        emboss = True
    )
    row_frame.prop(
        bpy.context.scene, 'mol_import_md_frame_step', 
        text = 'Step',
        emboss = True
    )
    row_frame.prop(
        bpy.context.scene, 'mol_import_md_frame_end', 
        text = 'End',
        emboss = True
    )
    col_main.prop(
        bpy.context.scene, 'mol_md_selection', 
        text = 'Import Filter', 
        emboss = True
    )
    col_main.separator()
    col_main.label(text="Custom Selections")
    row = col_main.row(align=True)
    
    row = row.split(factor = 0.9)
    row.template_list('MOL_UL_TrajectorySelectionListUI', 'A list', scene, 
                        "trajectory_selection_list", scene, "list_index", rows=3)
    col = row.column()
    col.operator('trajectory_selection_list.new_item', icon="ADD", text="")
    col.operator('trajectory_selection_list.delete_item', icon="REMOVE", text="")
    if scene.list_index >= 0 and scene.trajectory_selection_list:
        item = scene.trajectory_selection_list[scene.list_index]
        
        col = col_main.column(align=False)
        col.separator()
        
        col.prop(item, "name")
        col.prop(item, "selection")
