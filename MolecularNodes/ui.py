import bpy
from .pref import *
from .tools import property_exists
from . import nodes
from . import pkg
from . import load
from . import md
from . import assembly
import os



# operator that calls the function to import the structure from the PDB
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
        
        mol_object = load.molecule_rcsb(
            pdb_code=pdb_code,
            center_molecule=bpy.context.scene.mol_import_center, 
            del_solvent=bpy.context.scene.mol_import_del_solvent,
            include_bonds=bpy.context.scene.mol_import_include_bonds,
            starting_style=bpy.context.scene.mol_import_default_style
        )
        
        bpy.context.view_layer.objects.active = mol_object
        self.report({'INFO'}, message=f"Imported '{pdb_code}' as {mol_object.name}")
        
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
        include_bonds = bpy.context.scene.mol_import_include_bonds
        
        
        mol_object = load.molecule_local(
            file_path=file_path, 
            mol_name=bpy.context.scene.mol_import_local_name,
            include_bonds=include_bonds, 
            center_molecule=bpy.context.scene.mol_import_center, 
            del_solvent=bpy.context.scene.mol_import_del_solvent, 
            default_style=bpy.context.scene.mol_import_default_style, 
            setup_nodes=True
            )
        
        # return the good news!
        bpy.context.view_layer.objects.active = mol_object
        self.report({'INFO'}, message=f"Imported '{file_path}' as {mol_object.name}")
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
        file_top = bpy.context.scene.mol_import_md_topology
        file_traj = bpy.context.scene.mol_import_md_trajectory
        name = bpy.context.scene.mol_import_md_name
        
        mol_object, coll_frames = md.load_trajectory(
            file_top = file_top, 
            file_traj = file_traj, 
            name = name
        )
        n_frames = len(coll_frames.objects)
        
        nodes.create_starting_node_tree(
            obj = mol_object, 
            coll_frames = coll_frames, 
            starting_style = bpy.context.scene.mol_import_default_style
            )
        bpy.context.view_layer.objects.active = mol_object
        self.report({'INFO'}, message=f"Imported '{file_top}' as {mol_object.name} with {str(n_frames)} frames from '{file_traj}'.")
        
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
    row_name = col_main.row(align = True)
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
    row_import = col_main.row(align = True)
    row_import.prop(
        bpy.context.scene, 'mol_import_md_name', 
        text = "Name", 
        emboss = True
    )
    row_import.operator('mol.import_protein_md', text = "Load", icon_value = 30, emboss = True)
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

class MOL_OT_Default_Style(bpy.types.Operator):
    bl_idname = "mol.default_style"
    bl_label = "Change the default style."
    bl_description = "Change the default style of molecules on import."
    bl_options = {"REGISTER", "UNDO"}
    panel_display: bpy.props.IntProperty(name='panel_display', default = 0)

    @classmethod
    def poll(cls, context):
        return True

    def execute(self, context):
        bpy.context.scene.mol_import_default_style = self.panel_display
        return {"FINISHED"}


def default_style(layout, label, panel_display):
    op = layout.operator(
        'mol.default_style', 
        text = label, 
        emboss = True, 
        depress = (panel_display == bpy.context.scene.mol_import_default_style)
        )
    op.panel_display = panel_display

class MOL_MT_Default_Style(bpy.types.Menu):
    bl_label = ""
    bl_idname = "MOL_MT_Default_Style"
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout.column_flow(columns = 1)
        default_style(layout, 'Atoms', 0)
        default_style(layout, 'Ribbon', 1)
        default_style(layout, 'Ball and Stick', 2)

def MOL_PT_panel_ui(layout_function, ): 
    layout_function.label(text = "Import Options", icon = "MODIFIER")
    if not pkg.available():
        layout_function.operator('mol.install_dependencies', text = 'Install Packages')
    else:
        box = layout_function.box()
        grid = box.grid_flow(columns = 2)
        
        grid.prop(bpy.context.scene, 'mol_import_center', text = 'Centre Structre', icon_value=0, emboss=True)
        grid.prop(bpy.context.scene, 'mol_import_del_solvent', text = 'Delete Solvent', icon_value=0, emboss=True)
        grid.prop(bpy.context.scene, 'mol_import_include_bonds', text = 'Import Bonds', icon_value=0, emboss=True)
        grid.menu(
            'MOL_MT_Default_Style', 
            text = ['Atoms', 'Ribbon', 'Ball and Stick'][bpy.context.scene.mol_import_default_style])
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
        
        MOL_PT_panel_ui(self.layout, )


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

class MOL_OT_Assembly_Bio(bpy.types.Operator):
    bl_idname = "mol.assembly_bio"
    bl_label = "Build"
    bl_description = "Add Node to Build Biological Assembly"
    bl_options = {"REGISTER", "UNDO"}
    
    @classmethod
    def poll(cls, context):
        return True

    def execute(self, context):
        obj = context.active_object
        try:
            node_bio_assembly = assembly.create_biological_assembly_node(
                name = obj.name, 
                transform_dict = assembly.get_transformations_mmtf(obj['bio_transform_dict'])
            )
        except:
            node_bio_assembly = None
            self.report({'WARNING'}, message = 'Unable to detect biological assembly information.')
        
        if node_bio_assembly:
            mol_add_node(node_bio_assembly.name)
        
        return {"FINISHED"}


def menu_item_surface_custom(layout_function, label):
    op = layout_function.operator('mol.style_surface_custom', 
                                  text = label, 
                                  emboss = True, 
                                  depress = True)

def menu_item_color_chains(layout_function, label):
    op = layout_function.operator('mol.color_chains', 
                                  text = label, 
                                  emboss = True, 
                                  depress = True)

class MOL_OT_Color_Chain(bpy.types.Operator):
    bl_idname = "mol.color_chains"
    bl_label = "My Class Name"
    bl_description = "Create a custom node for coloring each chain individually."
    bl_options = {"REGISTER", "UNDO"}
    
    @classmethod
    def poll(cls, context):
        return True

    def execute(self, context):
        obj = context.active_object
        try:
            node_color_chain = nodes.chain_color(
                node_name = f"MOL_color_chains_{obj.name}", 
                input_list = obj['chain_id_unique']
            )
            mol_add_node(node_color_chain.name)
        except:
            self.report({'WARNING'}, message = 'Unable to detect chain information.')
        
        return {"FINISHED"}

def menu_chain_selection_custom(layout_function):
    obj = bpy.context.view_layer.objects.active
    label = 'Chain ' + str(obj.name)
    op = layout_function.operator(
        'mol.chain_selection_custom', 
        text = label, 
        emboss = True, 
        depress = True
    )

class MOL_OT_Chain_Selection_Custom(bpy.types.Operator):
    bl_idname = "mol.chain_selection_custom"
    bl_label = "Chain Selection"
    bl_description = "Add a custom node for selection all of the chains for this moledcule."
    bl_options = {"REGISTER", "UNDO"}
    
    @classmethod
    def poll(cls, context):
        return True
    
    def execute(self, context):
        obj = bpy.context.view_layer.objects.active
        node_chains = nodes.chain_selection(
            node_name = 'MOL_sel_' + str(obj.name) + "_chains", 
            input_list = obj['chain_id_unique']
            )
        
        mol_add_node(node_chains.name)
        
        return {"FINISHED"}


class MOL_MT_Add_Node_Menu_Properties(bpy.types.Menu):
    bl_idname = 'MOL_MT_ADD_NODE_MENU_PROPERTIES'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        # currently nothing for this menu in the panel

class MOL_MT_Add_Node_Menu_Color(bpy.types.Menu):
    bl_idname = 'MOL_MT_ADD_NODE_MENU_COLOR'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Set Color', 'MOL_color_set')
        layout.separator()
        menu_item_interface(layout, 'Goodsell Colors', 'MOL_color_goodsell')
        layout.separator()
        menu_item_interface(layout, 'Color by Atomic Number', 'MOL_color_atomic_number')
        menu_item_interface(layout, 'Color by Element', 'MOL_color_element')
        menu_item_color_chains(layout, 'Color by Chains')
        menu_item_interface(layout, 'Color Atomic', 'MOL_style_color')

class MOL_MT_Add_Node_Menu_Bonds(bpy.types.Menu):
    bl_idname = 'MOL_MT_ADD_NODE_MENU_BONDS'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Find Bonds', 'MOL_bonds_find')
        menu_item_interface(layout, 'Break Bonds', 'MOL_bonds_break')
        menu_item_interface(layout, 'Find Bonded Atoms', 'MOL_bonds_find_bonded')

class MOL_MT_Add_Node_Menu_Styling(bpy.types.Menu):
    bl_idname = 'MOL_MT_ADD_NODE_MENU_SYLING'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Atoms Cycles', 'MOL_style_atoms')
        menu_item_interface(layout, 'Atoms EEVEE', 'MOL_style_atoms_eevee')
        menu_item_interface(layout, 'Ribbon', 'MOL_style_ribbon')
        menu_item_interface(layout, 'Surface', 'MOL_style_surface_single')
        menu_item_surface_custom(layout, 'Surface Split Chains')
        menu_item_interface(layout, 'Ball and Stick', 'MOL_style_ball_and_stick')


class MOL_MT_Add_Node_Menu_Selections(bpy.types.Menu):
    bl_idname = 'MOL_MT_ADD_NODE_MENU_SELECTIONS'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Select Atoms', 'MOL_sel_atoms')
        menu_item_interface(layout, 'Separate Polymers', 'MOL_sel_sep_polymers')
        layout.separator()
        menu_chain_selection_custom(layout)
        layout.separator()
        menu_item_interface(layout, 'Atom Properties', 'MOL_sel_atom_propeties')
        menu_item_interface(layout, 'Atomic Number', 'MOL_sel_atomic_number')
        menu_item_interface(layout, 'Element Name', 'MOL_sel_element_name')
        layout.separator()
        menu_item_interface(layout, 'Res Properties', 'MOL_sel_res_properties')
        menu_item_interface(layout, 'Res Name', 'MOL_sel_res_name')
        menu_item_interface(layout, 'Res Name Nucleic', 'MOL_sel_res_name_nucleic')
        menu_item_interface(layout, 'Res ID', 'MOL_sel_res_id')
        menu_item_interface(layout, 'Res ID Range', 'MOL_sel_res_id_range')

class MOL_MT_Add_Node_Menu_Assembly(bpy.types.Menu):
    bl_idname = 'MOL_MT_ADD_NODE_MENU_ASSEMBLY'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        layout.operator("mol.assembly_bio", text = "Biological Assembly", emboss = True, depress=True)
        menu_item_interface(layout, 'Center Assembly', 'MOL_assembly_center')

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

class MOL_MT_Add_Node_Menu_DNA(bpy.types.Menu):
    bl_idname = 'MOL_MT_ADD_NODE_MENU_DNA'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Double Helix', 'MOL_dna_double_helix')
        menu_item_interface(layout, 'Bases', 'MOL_dna_bases')
        layout.separator()
        menu_item_interface(layout, 'Style Atoms Cyeles', 'MOL_dna_style_atoms')
        menu_item_interface(layout, 'Style Atoms EEVEE', 'MOL_dna_style_atoms_eevee')
        menu_item_interface(layout, 'Style Surface', 'MOL_dna_style_surface')
        menu_item_interface(layout, 'Style Ball and Stick', 'MOL_dna_style_ball_and_stick')

class MOL_MT_Add_Node_Menu_Animation(bpy.types.Menu):
    bl_idname = 'MOL_MT_ADD_NODE_MENU_ANIMATION'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Animate Frames', 'MOL_animate_frames')
        menu_item_interface(layout, 'Animate Value', 'MOL_animate_value')
        layout.separator()
        menu_item_interface(layout, 'Noise Position', 'MOL_noise_position')
        menu_item_interface(layout, 'Noise Field', 'MOL_noise_field')
        menu_item_interface(layout, 'Noise Repeat', 'MOL_noise_repeat')

class MOL_MT_Add_Node_Menu_Utilities(bpy.types.Menu):
    bl_idname = 'MOL_MT_ADD_NODE_MENU_UTILITIES'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Booelean Chain', 'MOL_utils_bool_chain')
        menu_item_interface(layout, 'Rotation Matrix', 'MOL_utils_rotation_matrix')
        menu_item_interface(layout, 'Curve Resample', 'MOL_utils_curve_resample')

class MOL_MT_Add_Node_Menu(bpy.types.Menu):
    bl_idname = "MOL_MT_ADD_NODE_MENU"
    bl_label = "Menu for Adding Nodes in GN Tree"

    @classmethod
    def poll(cls, context):
        return not (False)

    def draw(self, context):
        layout = self.layout.column_flow(columns=1)
        layout.operator_context = "INVOKE_DEFAULT"
        # layout.menu('MOL_MT_ADD_NODE_MENU_PROPERTIES', text='Properties', icon_value=201)
        layout.menu('MOL_MT_ADD_NODE_MENU_SYLING', text='Style', icon_value=77)
        layout.menu('MOL_MT_ADD_NODE_MENU_COLOR', text='Color', icon = 'COLORSET_07_VEC')
        layout.menu('MOL_MT_ADD_NODE_MENU_BONDS', text='Bonds', icon = 'FIXED_SIZE')
        layout.menu('MOL_MT_ADD_NODE_MENU_SELECTIONS', text='Selection', icon_value=256)
        layout.menu('MOL_MT_ADD_NODE_MENU_ANIMATION', text='Animation', icon_value=409)
        layout.menu('MOL_MT_ADD_NODE_MENU_ASSEMBLY', text='Assemblies', icon = 'GROUP_VERTEX')
        # layout.menu('MOL_MT_ADD_NODE_MENU_MEMBRANES', text='Membranes', icon_value=248)
        layout.menu('MOL_MT_ADD_NODE_MENU_DNA', text='DNA', icon='GP_SELECT_BETWEEN_STROKES')
        layout.menu('MOL_MT_ADD_NODE_MENU_UTILITIES', text='Utilities', icon_value=92)

def mol_add_node_menu(self, context):
    if ('GeometryNodeTree' == bpy.context.area.spaces[0].tree_type):
        layout = self.layout
        layout.menu('MOL_MT_ADD_NODE_MENU', text='Molecular Nodes', icon_value=88)