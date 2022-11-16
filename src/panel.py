import bpy
from . import open
from ..preferences import *
import os
import biotite.structure as struc

# operator that calls the function to import the structure frin tge PDB
class MOL_OT_Import_Protein_RCSB(bpy.types.Operator):
    bl_idname = "mol.import_protein_rcsb"
    bl_label = "import_protein_fetch_pdb"
    bl_description = "Download and open a structure from the Protein Data Bank"
    bl_options = {"REGISTER", "UNDO"}

    @classmethod
    def poll(cls, context):
        return not False

    def execute(self, context):
        mol = open.open_structure_rcsb(pdb_code = bpy.context.scene.mol_pdb_code)
        open.MOL_import_mol(
            mol = mol,
            mol_name = bpy.context.scene.mol_pdb_code,
            center_molecule = bpy.context.scene.mol_import_center,
            del_solvent = bpy.context.scene.mol_import_del_solvent, 
            include_bonds = bpy.context.scene.mol_import_include_bonds
            )
        return {"FINISHED"}

    def invoke(self, context, event):
        return self.execute(context)

# operator that calls the function to import the structure from a local file
class MOL_OT_Import_Protein_Local(bpy.types.Operator):
    bl_idname = "mol.import_protein_local"
    bl_label = "import_protein_local"
    bl_description = "Open a local structure file"
    bl_options = {"REGISTER", "UNDO"}

    @classmethod
    def poll(cls, context):
        return not False

    def execute(self, context):
        file_path = bpy.context.scene.mol_import_local_path
        file_ext = os.path.splitext(file_path)[1]
        include_bonds = bpy.context.scene.mol_import_include_bonds
        
        if file_ext == '.pdb':
            mol = open.open_structure_local_pdb(file_path, include_bonds)
        elif file_ext == '.pdbx' or file_ext == '.cif':
            mol = open.open_structure_local_pdbx(file_path, include_bonds)
            
        open.MOL_import_mol(
            mol = mol,
            mol_name = bpy.context.scene.mol_import_local_name,
            center_molecule = bpy.context.scene.mol_import_center,
            del_solvent = bpy.context.scene.mol_import_del_solvent, 
            include_bonds = include_bonds
            )
        return {"FINISHED"}

    def invoke(self, context, event):
        return self.execute(context)

def MOL_PT_panel_rcsb(layout_function, ):
    col_main = layout_function.column(heading = '', align = False)
    col_main.alert = False
    col_main.enabled = True
    col_main.active = True
    col_main.use_property_split = False
    col_main.use_property_decorate = False
    col_main.scale_x = 1.0
    col_main.scale_y = 1.0
    col_main.alignment = 'Expand'.upper()
    col_main.label(text = "Download from PDB", icon_value = 3)
    row_import = col_main.row()
    row_import.prop(bpy.context.scene, 'mol_pdb_code', text='PDB ID', icon_value=0, emboss=True)
    row_import.operator('mol.import_protein_rcsb', text='Download', icon_value=169, emboss=True, depress=False)

def MOL_PT_panel_local(layout_function, ):
    col_main = layout_function.column(heading = '', align = False)
    col_main.alert = False
    col_main.enabled = True
    col_main.active = True
    col_main.label(text = "Local Imported")
    row_name = col_main.row()
    row_name.prop(bpy.context.scene, 'mol_import_local_name', text = "Name", icon_value = 0, emboss = True)
    row_name.operator('mol.import_protein_local', text = "Load", icon_value = 30, emboss = True)
    row_import = col_main.row()
    row_import.prop(
        bpy.context.scene, 'mol_import_local_path', 
        text = "File path", 
        icon_value = 0, 
        emboss = True
    )

def MOL_PT_panel_md_traj(layout_function, ):
    col_main = layout_function.column(heading = '', align = False)
    col_main.alert = False
    col_main.enabled = True
    col_main.active = True
    col_main.label(text = "Import Molecular Dynamics Trajectories")
    row_topology = col_main.row(align = True)
    row_topology.prop(
        bpy.context.scene, 'mol_import_md_topology', 
        text = 'Topology', 
        icon_value = 458, 
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

class MOL_OT_Import_Method_Selection(bpy.types.Operator):
    bl_idname = "mol.import_method_selection"
    bl_label = "import_method"
    bl_description = "Change Structure Import Method"
    bl_options = {"REGISTER", "UNDO"}
    mol_interface_value: bpy.props.IntProperty(name = 'interface_value', description = '', default = 0, subtype = 'NONE')

    @classmethod
    def poll(cls, context):
        return not False

    def execute(self, context):
        bpy.context.scene.mol_import_panel_selection = self.mol_interface_value
        return {"FINISHED"}

    def invoke(self, context, event):
        return self.execute(context)

def MOL_change_import_interface(layout_function, label, interface_value, icon):
    op = layout_function.operator(
        'mol.import_method_selection', 
        text = label, 
        icon_value = icon, 
        emboss = True, 
        depress = interface_value == bpy.context.scene.mol_import_panel_selection
    )
    op.mol_interface_value = interface_value


def MOL_PT_panel_ui(layout_function, ): 
    layout_function.label(text = "Import Options", icon = "MODIFIER")
    box = layout_function.box()
    grid = box.grid_flow(columns = 2)
    
    grid.prop(bpy.context.scene, 'mol_import_center', text = 'Centre Structre', icon_value=0, emboss=True)
    grid.prop(bpy.context.scene, 'mol_import_del_solvent', text = 'Delete Solvent', icon_value=0, emboss=True)
    grid.prop(bpy.context.scene, 'mol_import_include_bonds', text = 'Import Bonds', icon_value=0, emboss=True)
    grid.label(text = "Default Style: Atoms")
    box = layout_function
    row = box.row(heading = '', align=True)
    row.alignment = 'EXPAND'
    row.enabled = True
    row.alert = False
    MOL_change_import_interface(row, 'PDB',           0,  72)
    MOL_change_import_interface(row, 'Local File',    1, 108)
    MOL_change_import_interface(row, 'MD Trajectory', 2, 487)
    
    layout_function = box.box()
    if bpy.context.scene.mol_import_panel_selection == 0:
        MOL_PT_panel_rcsb(layout_function)
    elif bpy.context.scene.mol_import_panel_selection == 1:
        MOL_PT_panel_local(layout_function)
    else:
        MOL_PT_panel_md_traj(layout_function)

class MOL_PT_panel(bpy.types.Panel):
    bl_label = 'Molecular Nodes'
    bl_idname = 'MOL_PT_panel'
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = 'scene'
    bl_order = 0
    bl_options = {'HEADER_LAYOUT_EXPAND'}
    bl_ui_units_x=0

    @classmethod
    def poll(cls, context):
        return not (False)

    def draw_header(self, context):
        layout = self.layout

    def draw(self, context):
        layout = self.layout
        layout_function = layout
        MOL_PT_panel_ui(layout_function, )
