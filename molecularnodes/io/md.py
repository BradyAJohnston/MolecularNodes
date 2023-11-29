"""
Importing molecular dynamics trajectories and associated files.
"""

__name__ = "MolecularNodes.trajectory"
__author__ = "Brady Johnston"

import bpy

try:
    import MDAnalysis as mda
except ImportError:
    HAS_mda = False
    import types

    class MockUniverse:
        pass

    mda = types.ModuleType("MDAnalysis")
    mda.Universe = MockUniverse

else:
    HAS_mda = True

from .mda import MDAnalysisSession
from .. import pkg

bpy.types.Scene.MN_import_md_topology = bpy.props.StringProperty(
    name='Topology',
    description='File path for the toplogy file for the trajectory',
    subtype='FILE_PATH',
    maxlen=0
)
bpy.types.Scene.MN_import_md_trajectory = bpy.props.StringProperty(
    name='Trajectory',
    description='File path for the trajectory file for the trajectory',
    subtype='FILE_PATH',
    maxlen=0
)
bpy.types.Scene.MN_import_md_name = bpy.props.StringProperty(
    name='Name',
    description='Name of the molecule on import',
    default='NewTrajectory',
    maxlen=0
)
bpy.types.Scene.MN_import_md_frame_start = bpy.props.IntProperty(
    name="Start",
    description="Frame start for importing MD trajectory",
    default=0
)
bpy.types.Scene.MN_import_md_frame_step = bpy.props.IntProperty(
    name="Step",
    description="Frame step for importing MD trajectory",
    default=1
)
bpy.types.Scene.MN_import_md_frame_stop = bpy.props.IntProperty(
    name="Stop",
    description="Frame stop for importing MD trajectory",
    default=499
)
bpy.types.Scene.MN_md_selection = bpy.props.StringProperty(
    name='Import Filter',
    description='Custom MDAnalysis selection string, removing unselecte atoms. See: "https://docs.mdanalysis.org/stable/documentation_pages/selections.html"',
    default='all'
)
bpy.types.Scene.MN_md_in_memory = bpy.props.BoolProperty(
    name='Memory',
    description='True will load all of the requested frames into the scene and into memory. False will stream the trajectory from a live MDAnalysis session',
    default=False
)
bpy.types.Scene.list_index = bpy.props.IntProperty(
    name="Index for trajectory selection list.",
    default=0
)

def load(
    top, 
    traj, 
    name = "NewTrajectory", 
    style = "spheres",
    selection:str = 'all', 
    start:int = 0, 
    step:int = 1,
    stop:int = 499, 
    subframes:int = 0, 
    custom_selections:dict = {}, 
    in_memory:bool = False
):
    universe = mda.Universe(top, traj)
    
    if in_memory:
        universe.transfer_to_memory(start=start, step=step, stop=stop)
        
    session = MDAnalysisSession()
    
    extra_selections = {}
    for sel in custom_selections:
        extra_selections[sel.name] = sel.selection
    mol = session.show(atoms=universe,
                            name=name,
                            style=style,
                            selection=selection,
                            custom_selections=extra_selections,
                            in_memory=in_memory
                            )
    
    return mol, universe


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
        if not pkg.is_current('MDAnalysis'):
            self.report({'ERROR'},
                        message="MDAnalysis is not installed. "
                                "Please install it to use this feature.")
            return {'CANCELLED'}
        top = scene.MN_import_md_topology
        traj = scene.MN_import_md_trajectory
        name = scene.MN_import_md_name

        mol, universe = load(
            top=top,
            traj = traj, 
            name = name,
            style = scene.MN_import_style,
            selection = scene.MN_md_selection, 
            start = scene.MN_import_md_frame_start, 
            stop = scene.MN_import_md_frame_stop, 
            step = scene.MN_import_md_frame_step, 
            custom_selections=scene.trajectory_selection_list,
            in_memory=scene.MN_md_in_memory
            
        )

        bpy.context.view_layer.objects.active = mol

        self.report(
            {'INFO'},
            message=f"Imported '{top}' as {name} "
                    f"with {str(universe.trajectory.n_frames)} "
                    f"frames from '{traj}'."
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

class MN_UL_TrajectorySelectionListUI(bpy.types.UIList):
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

def custom_selections(layout, scene):
    layout.label(text="Custom Selections")
    row = layout.row(align=True)
    
    row = row.split(factor = 0.9)
    row.template_list('MN_UL_TrajectorySelectionListUI', 'A list', scene, 
                        "trajectory_selection_list", scene, "list_index", rows=3)
    col = row.column()
    col.operator('trajectory_selection_list.new_item', icon="ADD", text="")
    col.operator('trajectory_selection_list.delete_item', icon="REMOVE", text="")
    if scene.list_index >= 0 and scene.trajectory_selection_list:
        item = scene.trajectory_selection_list[scene.list_index]
        
        col = layout.column(align=False)
        col.separator()
        
        col.prop(item, "name")
        col.prop(item, "selection")

def panel(layout, scene):
    layout.label(text = "Load MD Trajectories", icon='FILE_TICK')
    layout.separator()
    col = layout.column(align=True)
    row_import = col.row()
    row_import.prop(scene, 'MN_import_md_name')
    row_import.operator('mn.import_protein_md', text = "Load")
    col.separator()
    col.prop(scene, 'MN_import_md_topology')
    col.prop(scene, 'MN_import_md_trajectory')
    
    layout.separator()
    layout.label(text = "Options", icon = "MODIFIER")
    layout.prop(scene, "MN_import_style")
    layout.prop(scene, 'MN_md_selection')
    row_frame = layout.row(heading = "Frames", align = True)
    row_frame.prop(scene, 'MN_md_in_memory')
    row = row_frame.row(align=True)
    row.prop(scene, 'MN_import_md_frame_start')
    row.prop(scene, 'MN_import_md_frame_step')
    row.prop(scene, 'MN_import_md_frame_stop')
    row.enabled = scene.MN_md_in_memory
    custom_selections(layout, scene)