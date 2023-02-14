import bpy
from .pref import *
from .tools import property_exists
from . import nodes
from . import pkg
from . import load
from . import md
from . import assembly
import os,pathlib



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
        
        mol_object = load.molecule_local(
            file_path=file_path, 
            mol_name=bpy.context.scene.mol_import_local_name,
            include_bonds=bpy.context.scene.mol_import_include_bonds, 
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
        selection = bpy.context.scene.mol_md_selection
        md_start = bpy.context.scene.mol_import_md_frame_start
        md_step =  bpy.context.scene.mol_import_md_frame_step
        md_end =   bpy.context.scene.mol_import_md_frame_end
        del_solvent = bpy.context.scene.mol_import_del_solvent
        include_bonds = bpy.context.scene.mol_import_include_bonds
        
        mol_object, coll_frames = md.load_trajectory(
            file_top    = file_top, 
            file_traj   = file_traj, 
            md_start    = md_start,
            md_end      = md_end,
            md_step     = md_step,
            name        = name, 
            del_solvent = del_solvent, 
            selection   = selection,
            include_bonds=include_bonds
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
    col_main.label(text = "Download from PDB")
    row_import = col_main.row()
    row_import.prop(bpy.context.scene, 'mol_pdb_code', text='PDB ID')
    row_import.operator('mol.import_protein_rcsb', text='Download', icon='IMPORT')

def MOL_PT_panel_local(layout_function, ):
    col_main = layout_function.column(heading = '', align = False)
    col_main.alert = False
    col_main.enabled = True
    col_main.active = True
    col_main.label(text = "Open Local File")
    row_name = col_main.row(align = False)
    row_name.prop(bpy.context.scene, 'mol_import_local_name', text = "Name", icon_value = 0, emboss = True)
    row_name.operator('mol.import_protein_local', text = "Load", icon='FILE_TICK', emboss = True)
    row_import = col_main.row()
    row_import.prop(
        bpy.context.scene, 'mol_import_local_path', 
        text = "File path", 
        icon_value = 0, 
        emboss = True
    )

def MOL_PT_panel_md_traj(layout_function, scene):
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
    row.template_list('TrajectorySelectionListUI', 'A list', scene, 
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

def MOL_PT_panel_ui(layout_function, scene): 
    layout_function.label(text = "Import Options", icon = "MODIFIER")
    if not pkg.available():
        
        col_main = layout_function.column(heading = '', align = False)
        col_main.alert = False
        col_main.enabled = True
        col_main.active = True
        col_main.use_property_split = False
        col_main.use_property_decorate = False
        col_main.scale_x = 1.0
        col_main.scale_y = 1.0
        
        col_main.alignment = 'Expand'.upper()
        col_main.label(text = "Set PyPI Mirror")
        row_import = col_main.row()
        row_import.prop(bpy.context.scene, 'pypi_mirror',text='PyPI')
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
            MOL_PT_panel_md_traj(layout_function, scene)



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
        
        MOL_PT_panel_ui(self.layout, bpy.context.scene)


def mol_add_node(node_name):
    prev_context = bpy.context.area.type
    bpy.context.area.type = 'NODE_EDITOR'
    # actually invoke the operator to add a node to the current node tree
    # use_transform=True ensures it appears where the user's mouse is and is currently being moved
    # so the user can place it where they wish
    bpy.ops.node.add_node('INVOKE_DEFAULT', type='GeometryNodeGroup', use_transform=True)
    bpy.context.area.type = prev_context
    bpy.context.active_node.node_tree = bpy.data.node_groups[node_name]
    bpy.context.active_node.width = 200.0
    # checks to see if the node as a 'Material' property, and if it does, set MOL_atomic_materic as that property
    if (property_exists("bpy.data.node_groups[bpy.context.active_object.modifiers.active.node_group.name].nodes[bpy.context.active_node.name].inputs['Material'].default_value", globals(), locals())):
        mat = nodes.mol_base_material()
        bpy.data.node_groups[bpy.context.active_object.modifiers.active.node_group.name].nodes[bpy.context.active_node.name].inputs['Material'].default_value = bpy.data.materials[mat.name]
    
class MOL_OT_Add_Custom_Node_Group(bpy.types.Operator):
    bl_idname = "mol.add_custom_node_group"
    bl_label = "Add Custom Node Group"
    # bl_description = "Add Molecular Nodes custom node group."
    bl_options = {"REGISTER", "UNDO"}
    node_name: bpy.props.StringProperty(
        name = 'node_name', 
        description = '', 
        default = '', 
        subtype = 'NONE', 
        maxlen = 0
    )
    node_description: bpy.props.StringProperty(
        name = "node_description", 
        description="", 
        default="Add MolecularNodes custom node group.", 
        subtype="NONE"
    )

    @classmethod
    def poll(cls, context):
        return True
    
    @classmethod
    def description(cls, context, properties):
        return properties.node_description
    
    def execute(self, context):
        try:
            nodes.mol_append_node(self.node_name)
            mol_add_node(self.node_name)
        except RuntimeError:
            self.report({'ERROR'}, message='Failed to add node. Ensure you are not in edit mode.')
        return {"FINISHED"}
    
    def invoke(self, context, event):
        return self.execute(context)


def menu_item_interface(layout_function, 
                        label, 
                        node_name, 
                        node_description='Add custom MolecularNodes node group.'):
    op = layout_function.operator('mol.add_custom_node_group', text = label, emboss = True, depress=False)
    op.node_name = node_name
    op.node_description = node_description


class MOL_OT_Style_Surface_Custom(bpy.types.Operator):
    bl_idname = "mol.style_surface_custom"
    bl_label = "My Class Name"
    bl_description = "Create a split surface representation.\nGenerates an isosurface based on atomic vdw_radii. Each chain has its own separate surface representation"
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
    bl_description = "**PDB Downloaded Structures Only**\nAdds node to build biological assembly based on symmetry operations that are extraced from the structure file. Currently this is only supported for structures that were downloaded from the PDB"
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

def menu_residues_selection_custom(layout_function):
    obj = bpy.context.view_layer.objects.active
    label = 'Res ID'
    op = layout_function.operator(
        'mol.residues_selection_custom', 
        text = label, 
        emboss = True, 
        depress = True
    )

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
    bl_description = "Create a custom node for coloring each chain of a structure individually.\nRequires chain information to be available from the structure"
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
    bl_description = "Create a selection based on the chains.\nThis node is built on a per-molecule basis, taking into account the chain_ids that were detected. If no chain information is available this node will not work"
    bl_options = {"REGISTER", "UNDO"}
    
    @classmethod
    def poll(cls, context):
        return True
    
    def execute(self, context):
        obj = bpy.context.view_layer.objects.active
        node_chains = nodes.chain_selection(
            node_name = 'MOL_sel_' + str(obj.name) + "_chains", 
            input_list = obj['chain_id_unique'], 
            starting_value = 0,
            attribute = 'chain_id', 
            label_prefix = "Chain "
            )
        
        mol_add_node(node_chains.name)
        
        return {"FINISHED"}


class MOL_OT_Residues_Selection_Custom(bpy.types.Operator):
    bl_idname = "mol.residues_selection_custom"
    bl_label = "Multiple Residue Selection"
    bl_description = "Create a selection based on the provided residue strings.\nThis node is built on a per-molecule basis, taking into account the residues that were input. "
    bl_options = {"REGISTER", "UNDO"}

    input_resid_string: bpy.props.StringProperty(
        name="Select residue IDs: ",
        description="Enter a string value.",
        default="19,94,1-16"
    )

    @classmethod
    def poll(cls, context):
        return True
    
    def execute(self, context):
        obj = bpy.context.view_layer.objects.active
        node_residues = nodes.resid_multiple_selection(
            node_name = 'MOL_sel_residues', 
            input_resid_string = self.input_resid_string, 
            )
    
        
        mol_add_node(node_residues.name)
        return {"FINISHED"}

    def invoke(self, context, event):
        return context.window_manager.invoke_props_dialog(self)



def menu_ligand_selection_custom(layout_function):
    obj = bpy.context.view_layer.objects.active
    label = 'Ligands ' + str(obj.name)
    op = layout_function.operator(
        'mol.ligand_selection_custom', 
        text = label, 
        emboss = True, 
        depress = True
    )

class MOL_OT_Ligand_Selection_Custom(bpy.types.Operator):
    bl_idname = "mol.ligand_selection_custom"
    bl_label = "Ligand Selection"
    bl_description = "Create a selection based on the ligands.\nThis node is built on a per-molecule basis, taking into account the chain_ids that were detected. If no chain information is available this node will not work"
    bl_options = {"REGISTER", "UNDO"}
    
    @classmethod
    def poll(cls, context):
        return True
    
    def execute(self, context):
        obj = bpy.context.view_layer.objects.active
        node_chains = nodes.chain_selection(
            node_name = 'MOL_sel_' + str(obj.name) + "_ligands", 
            input_list = obj['ligands'], 
            starting_value = 100, 
            attribute = 'res_name', 
            label_prefix = ""
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
        menu_item_interface(layout, 'Set Color', 'MOL_color_set', 
                            "Sets a new color for the selected atoms")
        layout.separator()
        menu_item_interface(layout, 'Goodsell Colors', 'MOL_color_goodsell', 
                            "Adjusts the given colors to copy the 'Goodsell Style'.\n" + 
                            "Darkens the non-carbon atoms and keeps the carbon atoms the same color. " +
                            "Highlights differences without being too visually busy")
        layout.separator()
        menu_item_interface(layout, 'Color by Atomic Number', 'MOL_color_atomic_number', 
                            "Creates a color based on atomic_number field")
        menu_item_interface(layout, 'Color by Element', 'MOL_color_element', 
                            "Choose a color for each of the first 20 elements")
        menu_item_color_chains(layout, 'Color by Chains')
        menu_item_interface(layout, 'Color Atomic', 'MOL_style_color', 
                            "Choose a color for the most common elements in PDB structures")

class MOL_MT_Add_Node_Menu_Bonds(bpy.types.Menu):
    bl_idname = 'MOL_MT_ADD_NODE_MENU_BONDS'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Find Bonds', 'MOL_bonds_find', 
                            "Finds bonds between atoms based on distance.\n" + 
                            "Based on the vdw_radii for each point, finds other points within a certain radius to create a bond to. " + 
                            "Does not preserve the index for the points. Does not detect bond type")
        menu_item_interface(layout, 'Break Bonds', 'MOL_bonds_break', 
                            "Will delete a bond between atoms that already exists based on a distance cutoff")
        menu_item_interface(layout, 'Find Bonded Atoms', 'MOL_bonds_find_bonded', 
                            "Based on an initial selection, finds atoms which are within a certain number of bonds away")

class MOL_MT_Add_Node_Menu_Styling(bpy.types.Menu):
    bl_idname = 'MOL_MT_ADD_NODE_MENU_SYLING'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Atoms Cycles', 'MOL_style_atoms_cycles', 
                            'A sphere atom representation, visible ONLY in Cycles. Based on point-cloud rendering')
        menu_item_interface(layout, 'Atoms EEVEE', 'MOL_style_atoms_eevee', 
                            'A sphere atom representation, visible in EEVEE and Cycles. Based on mesh instancing which slows down viewport performance')
        menu_item_interface(layout, 'Ribbon Protein', 'MOL_style_ribbon_protein', 
                            'Create a ribbon mesh based off of the alpha-carbons of the structure')
        menu_item_interface(layout, 'Ribbon Nucleic', 'MOL_style_ribbon_nucleic', 
                            'Create a ribbon mesh and instanced cylinders for nucleic acids.')
        menu_item_interface(layout, 'Surface', 'MOL_style_surface_single', 
                            'Create a single joined surface representation.\n' +
                            'Generates an isosurface based on atomic vdw_radii. All chains are part of the same surface. Use "Surface Split Chains" ' + 
                            'to have a single surface per chain')
        menu_item_surface_custom(layout, 'Surface Split Chains')
        menu_item_interface(layout, 'Ball and Stick', 'MOL_style_ball_and_stick', 
                            'A style node to create ball and stick representation.\n' +
                            'Icospheres are instanced on atoms and cylinders for bonds. Bonds can be detected if they are not present in the structure')


class MOL_MT_Add_Node_Menu_Selections(bpy.types.Menu):
    bl_idname = 'MOL_MT_ADD_NODE_MENU_SELECTIONS'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Select Atoms', 'MOL_sel_atoms', 
                            "Separate atoms based on a selection field.\n" +
                            "Takes atoms and splits them into the selected atoms the inverted atoms, based on a selection field")
        menu_item_interface(layout, 'Separate Polymers', 'MOL_sel_sep_polymers', 
                            "Separate the Geometry into the different polymers.\n" + 
                            "Outputs for protein, nucleic & sugars")
        layout.separator()
        menu_chain_selection_custom(layout)
        menu_ligand_selection_custom(layout)
        layout.separator()
        menu_item_interface(layout, 'Atom Properties', 'MOL_sel_atom_propeties', 
                            "Create a selection based on the properties of the atom.\n" + 
                            "Fields for is_alpha_carbon, is_backbone, is_peptide, is_nucleic, is_solvent and is_carb")
        menu_item_interface(layout, 'Atomic Number', 'MOL_sel_atomic_number', 
                            "Create a selection if input value equal to the atomic_number field.")
        menu_item_interface(layout, 'Element Name', 'MOL_sel_element_name', 
                            "Create a selection of particular elements by name. Only first 20 elements supported")
        layout.separator()
        menu_item_interface(layout, 'Distance', 'MOL_sel_distance', 
                            "Create a selection based on the distance to a selected object.\n" + 
                            "The cutoff is scaled based on the objects scale and the 'Scale Cutoff' value.")
        menu_item_interface(layout, 'Slice', 'MOL_sel_slice', 
                            "Create a selection that is a slice along one of the XYZ axes, based on the position of an object.")
        layout.separator()
        menu_residues_selection_custom(layout)                        
        menu_item_interface(layout, 'Res ID Single', 'MOL_sel_res_id', 
                            "Create a selection if res_id matches input field")
        menu_item_interface(layout, 'Res ID Range', 'MOL_sel_res_id_range', 
                            "Create a selection if the res_id is within the given thresholds")
        menu_item_interface(layout, 'Res Name Peptide', 'MOL_sel_res_name', 
                            "Create a selection of particular amino acids by name")
        menu_item_interface(layout, 'Res Name Nucleic', 'MOL_sel_res_name_nucleic', 
                            "Create a selection of particular nucleic acids by name")
        menu_item_interface(layout, 'Res Whole', 'MOL_sel_res_whole', 
                            "Expand the selection to every atom in a residue, if any of those atoms are in the initial selection")
        menu_item_interface(layout, 'Res Atoms', 'MOL_sel_res_atoms', 
                            "Create a selection based on the atoms of a residue.\n" +
                            "Selections for CA, backbone atoms (N, C, O), sidechain and backbone")

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
        menu_item_interface(layout, 'Center Assembly', 'MOL_assembly_center', 
                            "Center the structure on the world origin based on bounding box")

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
        menu_item_interface(layout, 'Double Helix', 'MOL_dna_double_helix', 
                            "Create a DNA double helix from an input curve.\n" + 
                            "Takes an input curve and instances for the bases, returns instances of the bases in a double helix formation")
        menu_item_interface(layout, 'Bases', 'MOL_dna_bases', 
                            "Provide the DNA bases as instances to be styled and passed onto the Double Helix node")
        layout.separator()
        menu_item_interface(layout, 'Style Atoms Cyeles', 'MOL_dna_style_atoms_cycles', 
                            "Style the DNA bases with spheres only visible in Cycles")
        menu_item_interface(layout, 'Style Atoms EEVEE', 'MOL_dna_style_atoms_eevee', 
                            "Style the DNA bases with spheres visible in Cycles and EEVEE")
        menu_item_interface(layout, 'Style Surface', 'MOL_dna_style_surface', 
                            "Style the DNA bases with surface representation")
        menu_item_interface(layout, 'Style Ball and Stick', 'MOL_dna_style_ball_and_stick', 
                            "Style the DNA bases with ball and stick representation")

class MOL_MT_Add_Node_Menu_Animation(bpy.types.Menu):
    bl_idname = 'MOL_MT_ADD_NODE_MENU_ANIMATION'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Animate Frames', 'MOL_animate_frames', 
                            "Interpolate between frames of a trajectory." + 
                            "Given a collection of frames for a trajectory, this node interpolates between them from start to finish based on the Animate field taking a value from 0 to 1. The positions of the Atoms are then moved based on this field")
        menu_item_interface(layout, 'Animate Field', 'MOL_animate_field')
        menu_item_interface(layout, 'Animate Value', 'MOL_animate_value', 
                            "Animates between given start and end values, based on the input start and end frame of the timeline. Clamped will limit the output to the 'To Min' and 'To Max', while unclamped will continue to interpolate past these values. 'Smoother Step' will ease in and out of these values, with default being linear interpolation")
        layout.separator()
        menu_item_interface(layout, 'Res Wiggle', "MOL_animate_res_wiggle", 
                            "Wiggles the side chains of amino acids based on b_factor, adding movement to a structure.")
        menu_item_interface(layout, 'Res to Curve', "MOL_animate_res_to_curve", 
                            "Takes atoms and maps them along a curve, as a single long peptide chain.")
        layout.separator()
        menu_item_interface(layout, 'Noise Position', 'MOL_noise_position', 
                            "Generate 3D noise field based on the position attribute")
        menu_item_interface(layout, 'Noise Field', 'MOL_noise_field', 
                            "Generate a 3D noise field based on the given field")
        menu_item_interface(layout, 'Noise Repeat', 'MOL_noise_repeat', 
                            "Generate a 3D noise field that repeats, based on the given field")

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
