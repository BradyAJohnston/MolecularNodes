import bpy
from . import open

# operator that calls the function to import the structure
class MOL_OT_Import_Protein_RCSB(bpy.types.Operator):
    bl_idname = "mol.import_protein_rcsb"
    bl_label = "import_protein_fetch_pdb"
    bl_description = "Download and open a structure from the Protein Data Bank"
    bl_options = {"REGISTER", "UNDO"}

    @classmethod
    def poll(cls, context):
        return not False

    def execute(self, context):
        open.import_protein_pdb(
            pdb_code = bpy.context.scene.mol_pdb_code, 
            center_molecule = bpy.context.scene.mol_import_center,
            del_solvent = bpy.context.scene.mol_import_del_solvent
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
    col_main.label(text = "Testing Again", icon_value = 3)
    row_import = col_main.row()
    row_import.prop(bpy.context.scene, 'mol_pdb_code', text='PDB ID', icon_value=0, emboss=True)
    row_import.operator('mol.import_protein_rcsb', text='Download', icon_value=169, emboss=True, depress=False)
    row_options = col_main.row()
    row_options.prop(bpy.context.scene, 'mol_import_center', text='Centre Structre', icon_value=0, emboss=True)
    row_options.prop(bpy.context.scene, 'mol_import_del_solvent', text='Delete Solvent', icon_value=0, emboss=True)

def MOL_PT_panel_local(layout_function, ):
    col_main = layout_function.column(heading = '', align = False)
    col_main.alert = False
    col_main.enabled = True
    col_main.active = True
    col_main.label(text = "Local Imported")
    row_import = col_main.row()
    row_import.prop(
        bpy.context.scene, 'mol_import_local_path', 
        text = "File path", 
        icon_value = 0, 
        emboss = True
    )

class MOL_OT_Import_Method_Selection(bpy.types.Operator):
    bl_idname = 'mol.import_method_selection', 
    bl_label = "import_method_selection", 
    bl_description = "Change Structure Import Method", 
    bl_options = {"REGISTER", "UNDO"}, 
    interface_value: bpy.props.IntProperty(name = 'interface_value', description = '', default = 0, subtype = 'NONE')

    @classmethod
    def poll(cls, context):
        return not False

    def execute(self, context):
        bpy.context.scene.mol_import_panel_selection = self.interface_value
        return {"FINISHED"}

    def invoke(self, context, event):
        return self.execute(context)

def MOL_change_import_interface(layout_function, label, interface_value, icon):
    op = layout_function.operator(
        'mol.import_method_selection', 
        text = label, 
        icon_value = icon, 
        emboss = True, 
        depress = (interface_value == bpy.context.scene.mol_import_panel_selection)
    )
    #op.mol_import_method_selection = interface_value

def MOL_PT_panel_ui(layout_function, ): 
    row = layout_function.row()
    MOL_change_import_interface(layout_function, 'PDB',           0,  72)
    MOL_change_import_interface(layout_function, 'Local File',    1, 108)
    MOL_change_import_interface(layout_function, 'MD Trajectory', 2, 487)

    if True: #bpy.context.scene.mol_import_panel_selection == 0:
        MOL_PT_panel_rcsb(layout_function)
    else: #bpy.context.scene.mol_import_panel_selection == 1:
        MOL_PT_panel_local(layout_function)
    #else:
    #    MOL_PT_panel_md(layout_function)

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