import bpy
from . import open
from ..preferences import *
from .tools import property_exists
from . import nodes
from . import md
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
        pdb_code = bpy.context.scene.mol_pdb_code
        mol = open.open_structure_rcsb(pdb_code = pdb_code)
        mol_object = open.MOL_import_mol(
            mol = mol,
            mol_name = pdb_code,
            center_molecule = bpy.context.scene.mol_import_center,
            del_solvent = bpy.context.scene.mol_import_del_solvent, 
            include_bonds = bpy.context.scene.mol_import_include_bonds
            )
        
        nodes.create_starting_node_tree(mol_object)
        self.report({'INFO'}, message='Successfully Imported '+ pdb_code + ' as ' + mol_object.name)
        
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
        file_path = os.path.abspath(file_path)
        file_ext = os.path.splitext(file_path)[1]
        include_bonds = bpy.context.scene.mol_import_include_bonds
        
        if file_ext == '.pdb':
            mol = open.open_structure_local_pdb(file_path, include_bonds)
        elif file_ext == '.pdbx' or file_ext == '.cif':
            mol = open.open_structure_local_pdbx(file_path, include_bonds)
        
        mol_name = bpy.context.scene.mol_import_local_name
        mol_object = open.MOL_import_mol(
            mol = mol,
            mol_name = mol_name,
            center_molecule = bpy.context.scene.mol_import_center,
            del_solvent = bpy.context.scene.mol_import_del_solvent, 
            include_bonds = include_bonds
            )
        # setup the required initial node tree on the object 
        nodes.create_starting_node_tree(mol_object)
        # return the good news!
        self.report({'INFO'}, message='Successfully Imported '+ file_path + " as " + mol_object.name)
        return {"FINISHED"}

    def invoke(self, context, event):
        return self.execute(context)

class MOL_OT_Import_Protein_MD(bpy.types.Operator):
    bl_idname = "mol.import_protein_md"
    bl_label = "Import Protein MD"
    bl_description = "Load molecular dynamics trajectory"
    bl_options = {"REGISTER", "UNDO"}

    @classmethod
    def poll(cls, context):
        return True

    def execute(self, context):
        mol_object = md.load_trajectory(
            file_top = bpy.context.scene.mol_import_md_topology, 
            file_traj = bpy.context.scene.mol_import_md_trajectory, 
            name = bpy.context.scene.mol_import_md_name
        )
        
        nodes.create_starting_node_tree(mol_object)
        self.report({'SUCCESS'}, message='Successfully Imported Trajectory')
        
        return {"FINISHED"}


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
    col_main.label(text = "Open Local File")
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
    row_import = col_main.row(align = True)
    row_import.prop(
        bpy.context.scene, 'mol_import_md_name', 
        text = "Name", 
        emboss = True
    )
    row_import.operator('mol.import_protein_md', text = "Load", icon_value = 30, emboss = True)
    

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


def mol_add_node(node_name):
    prev_context = bpy.context.area.type
    bpy.context.area.type = 'NODE_EDITOR'
    bpy.ops.node.add_node('INVOKE_DEFAULT', type='GeometryNodeGroup', use_transform=True)
    bpy.context.area.type = prev_context
    bpy.context.active_node.node_tree = bpy.data.node_groups[node_name]
    bpy.context.active_node.width = 200.0
    if (property_exists("bpy.data.node_groups[bpy.context.active_object.modifiers.active.node_group.name].nodes[bpy.context.active_node.name].inputs['Material'].default_value", globals(), locals())):
        mat = nodes.mol_base_material()
        bpy.data.node_groups[bpy.context.active_object.modifiers.active.node_group.name].nodes[bpy.context.active_node.name].inputs['Material'].default_value = bpy.data.materials[mat.name]
    
class MOL_OT_Add_Custom_Node_Group(bpy.types.Operator):
    bl_idname = "mol.add_custom_node_group"
    bl_label = "Add Custom Node Group"
    bl_description = "Add Molecular Nodes custom node group."
    bl_options = {"REGISTER", "UNDO"}
    node_name: bpy.props.StringProperty(
        name = 'node_name', 
        description = '', 
        default = '', 
        subtype = 'NONE', 
        maxlen = 0
    )

    @classmethod
    def poll(cls, context):
        return True

    def execute(self, context):
        try:
            nodes.mol_append_node(self.node_name)
            mol_add_node(self.node_name)
        except RuntimeError:
            self.report({'ERROR'}, message='Failed to add node. Ensure you are not in edit mode.')
        return {"FINISHED"}
    
    def invoke(self, context, event):
        return self.execute(context)


def menu_item_interface(layout_function, label, node):
    op = layout_function.operator('mol.add_custom_node_group', 
                                  text = label, 
                                  emboss = True, depress=False)
    op.node_name = node


class MOL_OT_Style_Surface_Custom(bpy.types.Operator):
    bl_idname = "mol.style_surface_custom"
    bl_label = "My Class Name"
    bl_description = "Create a surface representation for each chain."
    bl_options = {"REGISTER", "UNDO"}
    
    @classmethod
    def poll(cls, context):
        return True

    def execute(self, context):
        obj = context.active_object
        try:
            node_surface = nodes.create_custom_surface(
                name = 'MOL_style_surface_' + obj.name + '_split', 
                n_chains = len(obj['chain_id_unique'])
            )
        except:
            node_surface = nodes.mol_append_node('MOL_style_surface_single')
            self.report({'WARNING'}, message = 'Unable to detect number of chains.')
        mol_add_node(node_surface.name)
        
        return {"FINISHED"}


def menu_item_surface_custom(layout_function, label):
    op = layout_function.operator('mol.style_surface_custom', 
                                  text = label, 
                                  emboss = True, 
                                  depress = True)


class MOL_MT_Add_Node_Menu_Properties(bpy.types.Menu):
    bl_idname = 'MOL_MT_ADD_NODE_MENU_PROPERTIES'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Atomic Properties', 'MOL_prop_atomic')
        menu_item_interface(layout, 'AA  Atom Positions', 'MOL_prop_AA_pos')
        menu_item_interface(layout, 'AA  Unique Number', 'MOL_prop_unique_AA_num')
        menu_item_interface(layout, 'Resample Curve', 'MOL_utils_curve_resample')
        menu_item_interface(layout, 'Radii Rescale', 'MOL_prop_scale_radii')
        menu_item_interface(layout, 'Radii Lookup', 'MOL_prop_radii')
        menu_item_interface(layout, 'Find Bonds', 'MOL_prop_find_bonds')

class MOL_MT_Add_Node_Menu_Styling(bpy.types.Menu):
    bl_idname = 'MOL_MT_ADD_NODE_MENU_SYLING'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Ribbon', 'MOL_style_ribbon')
        menu_item_interface(layout, 'Surface', 'MOL_style_surface_single')
        menu_item_surface_custom(layout, 'Surface Split Chains')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')


class MOL_MT_Add_Node_Menu_Selections(bpy.types.Menu):
    bl_idname = 'MOL_MT_ADD_NODE_MENU_SELECTIONS'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')

class MOL_MT_Add_Node_Menu_Membranes(bpy.types.Menu):
    bl_idname = 'MOL_MT_ADD_NODE_MENU_MEMBRANES'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')

class MOL_MT_Add_Node_Menu_DNA(bpy.types.Menu):
    bl_idname = 'MOL_MT_ADD_NODE_MENU_DNA'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')

class MOL_MT_Add_Node_Menu_Animation(bpy.types.Menu):
    bl_idname = 'MOL_MT_ADD_NODE_MENU_ANIMATION'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')

class MOL_MT_Add_Node_Menu_Utilities(bpy.types.Menu):
    bl_idname = 'MOL_MT_ADD_NODE_MENU_UTILITIES'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')
        menu_item_interface(layout, 'Setup Atomic Properties', 'MOL_prop_setup')

class MOL_MT_Add_Node_Menu(bpy.types.Menu):
    bl_idname = "MOL_MT_ADD_NODE_MENU"
    bl_label = "Menu for Ading Nodes in GN Tree"

    @classmethod
    def poll(cls, context):
        return not (False)

    def draw(self, context):
        layout = self.layout.column_flow(columns=1)
        layout.operator_context = "INVOKE_DEFAULT"
        layout.menu('MOL_MT_ADD_NODE_MENU_PROPERTIES', text='Properties', icon_value=201)
        layout.menu('MOL_MT_ADD_NODE_MENU_SYLING', text='Styling', icon_value=77)
        layout.menu('MOL_MT_ADD_NODE_MENU_SELECTIONS', text='Selections', icon_value=256)
        layout.menu('MOL_MT_ADD_NODE_MENU_MEMBRANES', text='Membranes', icon_value=248)
        layout.menu('MOL_MT_ADD_NODE_MENU_DNA', text='DNA', icon_value=206)
        layout.menu('MOL_MT_ADD_NODE_MENU_ANIMATION', text='Animation', icon_value=409)
        layout.menu('MOL_MT_ADD_NODE_MENU_UTILITIES', text='Utilities', icon_value=92)

def mol_add_node_menu(self, context):
    if ('GeometryNodeTree' == bpy.context.area.spaces[0].tree_type):
        layout = self.layout
        layout.menu('MOL_MT_ADD_NODE_MENU', text='Molecular Nodes', icon_value=88)