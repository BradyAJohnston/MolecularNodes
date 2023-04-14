import bpy
import numpy as np
from . import data
from . import coll
from .load import create_object, add_attribute
import warnings
from .md import load_trajectory


class SoluteSelection(bpy.types.PropertyGroup):
    """Group of properties for custom selections for solute import."""
    
    name: bpy.props.StringProperty(
        name="Solute Name", 
        description="Attribute", 
        default="name_solute"
    )
    
    selection: bpy.props.StringProperty(
        name="Solute Selection String", 
        description="String that provides a selection through MDAnalysis", 
        default = "name CA"
    )

class SolventGroupSelection(bpy.types.PropertyGroup):
    """Group of properties for custom selections for solvent groups import."""
    
    name: bpy.props.StringProperty(
        name="Solvent Group Name", 
        description="Attribute", 
        default="name_solvent"
    )
    
    selection: bpy.props.StringProperty(
        name="Solvent Group Selection String", 
        description="String that provides a selection through MDAnalysis", 
        default = "name CA"
    )

    shell_number: bpy.props.IntProperty(
        name="Shell Number Input", 
        description="Shell Count ", 
        default = 0
    )


class MOL_UL_SoluteSelectionListUI(bpy.types.UIList):
    """UI List for Solute """
    
    def draw_item(self, context, layout, data, item, 
                  icon, active_data, active_propname, index):
        custom_icon = "VIS_SEL_11"
        
        if self.layout_type in {'DEFAULT', 'COMPACT'}:
            layout.label(text = item.name, icon = custom_icon)
        
        elif self.layout_type in {'GRID'}:
            layout.alignment = 'CENTER'
            layout.label(text = "", icon = custom_icon)

class MOL_UL_SolventGroupSelectionListUI(bpy.types.UIList):
    """UI List Solvent Group"""
    
    def draw_item(self, context, layout, data, item, 
                  icon, active_data, active_propname, index):
        custom_icon = "VIS_SEL_11"
        
        if self.layout_type in {'DEFAULT', 'COMPACT'}:
            layout.label(text = item.name, icon = custom_icon)
        
        elif self.layout_type in {'GRID'}:
            layout.alignment = 'CENTER'
            layout.label(text = "", icon = custom_icon)


class SoluteSelection_OT_NewItem(bpy.types.Operator):
    """Add a new custom selection to the list."""
    
    bl_idname = "solute_input.new_item"
    bl_label = "+"
    
    def execute(self, context):
        context.scene.solute_input.add()
        return {'FINISHED'}

class SoluteSelection_OT_DeleteIem(bpy.types.Operator):
    
    bl_idname = "solute_input.delete_item"
    bl_label = "-"
    
    @classmethod
    def poll(cls, context):
        return context.scene.solute_input
    def execute(self, context):
        my_list = context.scene.solute_input
        index = context.scene.solute_index
        
        my_list.remove(index)
        context.scene.solute_index = min(max(0, index - 1), len(my_list) - 1)
        
        return {'FINISHED'}
    

class SolventGroupSelection_OT_NewItem(bpy.types.Operator):
    """Add a new custom selection to the list."""
    
    bl_idname = "solvent_groups_list.new_item"
    bl_label = "+"
    
    def execute(self, context):
        context.scene.solvent_groups_list.add()
        return {'FINISHED'}

class SolventGroupSelection_OT_DeleteIem(bpy.types.Operator):
    
    bl_idname = "solvent_groups_list.delete_item"
    bl_label = "-"
    
    @classmethod
    def poll(cls, context):
        return context.scene.solvent_groups_list
    def execute(self, context):
        my_list = context.scene.solvent_groups_list
        index = context.scene.solvent_groups_list_index
        
        my_list.remove(index)
        context.scene.solvent_groups_list_index = min(max(0, index - 1), len(my_list) - 1)
        
        return {'FINISHED'}



def build_selections(

                    file_top, 
                    file_traj,
                    frame,
                    world_scale = 0.01, 
                    include_bonds = False, 
                    del_solvent = False,
                    solute=None, 
                    solvent=None,
                    name="Default"
                    
                    ):
    
    import MDAnalysis as mda

    from collections import namedtuple
    
    from solvation_analysis.solute import Solute

    # create mdanalysis universe
    # initially load in the trajectory
    if file_traj == "":
        univ = mda.Universe(file_top)
    else:
        univ = mda.Universe(file_top, file_traj)
        
    # separate the trajectory, separate to the topology or the subsequence selections
    traj = univ.trajectory[frame]

    # define solute atom group
    solute_selection = univ.select_atoms(solute[0].selection)
    # define solvent atom group
    solvent_selections = {}
    for solvent_group in solvent:
        solvent_selections[solvent_group.name] = univ.select_atoms(solvent_group.selection)

    # instantiate solution
    solution = Solute.from_atoms(solute_selection, solvent_selections)

    # run the solvation analysis.
    solution.run()

    # make dictionary for shell selection
    shell_selection = {}
    for solvent_group in solvent:
        shell_selection[solvent_group.name] = solvent_group.shell_number

    # Get shell info from solution dataframe
    shells = solution.speciation.get_shells(shell_selection)

    # shells is a dataframe with a column named "solute_ix" that we need

    frames =shells.index.get_level_values(0).to_list()
    centers = shells.index.get_level_values(1).to_list()
    frame_centers_dict=list(zip(frames,centers))


    # we will need something that mimics the custom selection object
    # it needs to have a name attribute and a selection attribute
    # we will use a Python named tuple
    mock_selection = namedtuple('mock_selection', ['name', 'selection'])

    shells = []
    for i in frame_centers_dict:
        if i[0]==frame:
            # save the AtomGroup
            shell_group = solution.get_shell(solute_index=i[1], frame=0).indices

            # Create string selections from shell_group
            new_selection_string = "index " + " ".join(f"{index}" for index in shell_group)

            # Create a custom selection object
            shells.append(mock_selection(f"shell_{i}", new_selection_string))


    # call load_trajectory
    mol_object,coll_frame=load_trajectory(file_top, file_traj, frame, frame+1, 1, world_scale, include_bonds, del_solvent=False, selection="", name=name, custom_selections=shells)

    return mol_object,coll_frame