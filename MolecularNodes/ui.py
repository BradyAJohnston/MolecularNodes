import bpy
from . import nodes
from . import pkg
from . import load
from . import md
from . import density
from . import star
from . import esmfold

# operator that calls the function to import the structure from the PDB
class MN_OT_Import_Protein_RCSB(bpy.types.Operator):
    bl_idname = "mn.import_protein_rcsb"
    bl_label = "import_protein_fetch_pdb"
    bl_description = "Download and open a structure from the Protein Data Bank"
    bl_options = {"REGISTER", "UNDO"}

    @classmethod
    def poll(cls, context):
        return not False

    def execute(self, context):
        pdb_code = bpy.context.scene.MN_pdb_code
        
        MN_object = load.molecule_rcsb(
            pdb_code=pdb_code,
            center_molecule=bpy.context.scene.MN_import_center, 
            del_solvent=bpy.context.scene.MN_import_del_solvent,
            include_bonds=bpy.context.scene.MN_import_include_bonds,
            starting_style=bpy.context.scene.MN_import_default_style,
            cache_dir=bpy.context.scene.MN_cache_dir
        )
        
        bpy.context.view_layer.objects.active = MN_object
        self.report({'INFO'}, message=f"Imported '{pdb_code}' as {MN_object.name}")
        
        return {"FINISHED"}

    def invoke(self, context, event):
        return self.execute(context)
    

# operator that calls the function to import the structure from a local file
class MN_OT_Import_Protein_Local(bpy.types.Operator):
    bl_idname = "mn.import_protein_local"
    bl_label = "import_protein_local"
    bl_description = "Open a local structure file"
    bl_options = {"REGISTER", "UNDO"}

    @classmethod
    def poll(cls, context):
        return not False

    def execute(self, context):
        file_path = bpy.context.scene.MN_import_local_path
        
        MN_object = load.molecule_local(
            file_path=file_path, 
            MN_name=bpy.context.scene.MN_import_local_name,
            include_bonds=bpy.context.scene.MN_import_include_bonds, 
            center_molecule=bpy.context.scene.MN_import_center, 
            del_solvent=bpy.context.scene.MN_import_del_solvent, 
            default_style=bpy.context.scene.MN_import_default_style, 
            setup_nodes=True
            )
        
        # return the good news!
        bpy.context.view_layer.objects.active = MN_object
        self.report({'INFO'}, message=f"Imported '{file_path}' as {MN_object.name}")
        return {"FINISHED"}

    def invoke(self, context, event):
        return self.execute(context)



def MN_PT_panel_rcsb(layout_function, ):
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
    col_main.prop(
        bpy.context.scene,
        'MN_cache_dir',
        text = 'Cache dir')
    row_import = col_main.row()
    row_import.prop(bpy.context.scene, 'MN_pdb_code', text='PDB ID')
    row_import.operator('mn.import_protein_rcsb', text='Download', icon='IMPORT')


def MN_PT_panel_local(layout_function, ):
    col_main = layout_function.column(heading = '', align = False)
    col_main.alert = False
    col_main.enabled = True
    col_main.active = True
    col_main.label(text = "Open Local File")
    row_name = col_main.row(align = False)
    row_name.prop(bpy.context.scene, 'MN_import_local_name', 
                    text = "Name", icon_value = 0, emboss = True)
    row_name.operator('mn.import_protein_local', text = "Load", 
                        icon='FILE_TICK', emboss = True)
    row_import = col_main.row()
    row_import.prop(
        bpy.context.scene, 'MN_import_local_path', 
        text = "File path", 
        icon_value = 0, 
        emboss = True
    )



class MN_OT_Import_Method_Selection(bpy.types.Operator):
    bl_idname = "mn.import_method_selection"
    bl_label = "import_method"
    bl_description = "Change Structure Import Method"
    bl_options = {"REGISTER", "UNDO"}
    MN_interface_value: bpy.props.IntProperty(
        name = 'interface_value', 
        description = '', 
        default = 0, 
        subtype = 'NONE'
        )

    @classmethod
    def poll(cls, context):
        return not False

    def execute(self, context):
        bpy.context.scene.MN_import_panel_selection = self.MN_interface_value
        return {"FINISHED"}

    def invoke(self, context, event):
        return self.execute(context)

def MN_change_import_interface(layout_function, label, interface_value, icon):
    if isinstance(icon, str):
        op = layout_function.operator(
            'mn.import_method_selection', 
            text = label, 
            icon = icon, 
            emboss = True, 
            depress = interface_value == bpy.context.scene.MN_import_panel_selection
        )
    elif isinstance(icon, int):
        op = layout_function.operator(
            'mn.import_method_selection', 
            text = label, 
            icon_value = icon, 
            emboss = True, 
            depress = interface_value == bpy.context.scene.MN_import_panel_selection
        )
    op.MN_interface_value = interface_value

class MN_OT_Default_Style(bpy.types.Operator):
    bl_idname = "mn.default_style"
    bl_label = "Change the default style."
    bl_description = "Change the default style of molecules on import."
    bl_options = {"REGISTER", "UNDO"}
    panel_display: bpy.props.IntProperty(name='panel_display', default = 0)

    @classmethod
    def poll(cls, context):
        return True

    def execute(self, context):
        bpy.context.scene.MN_import_default_style = self.panel_display
        return {"FINISHED"}

def default_style(layout, label, panel_display):
    op = layout.operator(
        'mn.default_style', 
        text = label, 
        emboss = True, 
        depress = (panel_display == bpy.context.scene.MN_import_default_style)
        )
    op.panel_display = panel_display

class MN_MT_Default_Style(bpy.types.Menu):
    bl_label = ""
    bl_idname = "MN_MT_Default_Style"
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout.column_flow(columns = 1)
        default_style(layout, 'Atoms', 0)
        default_style(layout, 'Cartoon', 1)
        default_style(layout, 'Ribbon', 2)
        default_style(layout, 'Ball and Stick', 3)

def MN_PT_panel_ui(layout_function, scene): 
    layout_function.label(text = "Import Options", icon = "MODIFIER")
    box = layout_function.box()
    grid = box.grid_flow(columns = 2)
    
    grid.prop(bpy.context.scene, 'MN_import_center', 
                text = 'Centre Structure', icon_value=0, emboss=True)
    grid.prop(bpy.context.scene, 'MN_import_del_solvent', 
                text = 'Delete Solvent', icon_value=0, emboss=True)
    grid.prop(bpy.context.scene, 'MN_import_include_bonds', 
                text = 'Import Bonds', icon_value=0, emboss=True)
    grid.menu(
        'MN_MT_Default_Style', 
        text = ['Atoms', 'Cartoon', 'Ribbon', 'Ball and Stick'][
            bpy.context.scene.MN_import_default_style
            ])
    panel = layout_function
    # row = panel.row(heading = '', align=True)
    row = panel.grid_flow(row_major = True, columns = 3, align = True)
    row.alignment = 'EXPAND'
    row.enabled = True
    row.alert = False
    
    
    MN_change_import_interface(row, 'PDB',           0,  "URL")
    MN_change_import_interface(row, 'ESMFold',       1,  "URL")
    MN_change_import_interface(row, 'Local File',    2, 108)
    MN_change_import_interface(row, 'MD Trajectory', 3, 487)
    MN_change_import_interface(row, 'EM Map', 4, 'LIGHTPROBE_CUBEMAP')
    MN_change_import_interface(row, 'Star File',     5, 487)
    
    panel_selection = bpy.context.scene.MN_import_panel_selection
    col = panel.column()
    box = col.box()
    
    if panel_selection == 0:
        row = layout_function.row()
        if not pkg.is_current('biotite'):
            box.enabled = False
            box.alert = True
            box.label(text = "Please install biotite in the addon preferences.")
        
        MN_PT_panel_rcsb(box)
    elif panel_selection == 1:
        if not pkg.is_current('biotite'):
            box.enabled = False
            box.alert = True
            box.label(text = "Please install biotite in the addon preferences.")
        esmfold.panel(box)
    elif panel_selection == 2:
        if not pkg.is_current('biotite'):
            box.enabled = False
            box.alert = True
            box.label(text = "Please install biotite in the addon preferences.")
        MN_PT_panel_local(box)
    elif panel_selection == 3:
        if not pkg.is_current('MDAnalysis'):
            box.enabled = False
            box.alert = True
            box.label(text = "Please install MDAnalysis in the addon preferences.")
            
        md.panel(box, scene)
    elif panel_selection == 4:
        if not pkg.is_current('mrcfile'):
            box.enabled = False
            box.alert = True
            box.label(text = "Please intall 'mrcfile' in the addon preferences.")
        density.panel(box, scene)
    elif panel_selection == 5:
        for name in ['starfile', 'eulerangles']:
            if not pkg.is_current(name):
                box.enabled = False
                box.alert = True
                box.label(text = f"Please install '{name}' in the addon preferences.")
        star.panel(box, scene)


class MN_PT_panel(bpy.types.Panel):
    bl_label = 'Molecular Nodes'
    bl_idname = 'MN_PT_panel'
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
        
        MN_PT_panel_ui(self.layout, bpy.context.scene)

def MN_add_node(node_name, label: str = '', show_options = False):
    prev_context = bpy.context.area.type
    bpy.context.area.type = 'NODE_EDITOR'
    # actually invoke the operator to add a node to the current node tree
    # use_transform=True ensures it appears where the user's mouse is and is currently 
    # being moved so the user can place it where they wish
    bpy.ops.node.add_node(
        'INVOKE_DEFAULT', 
        type='GeometryNodeGroup', 
        use_transform=True
        )
    bpy.context.area.type = prev_context
    node = bpy.context.active_node
    node.node_tree = bpy.data.node_groups[node_name]
    node.width = 200.0
    node.show_options = show_options
    if label != '':
        node.label = label
    
    # if added node has a 'Material' input, set it to the default MN material
    input_mat = bpy.context.active_node.inputs.get('Material')
    if input_mat:
        input_mat.default_value = nodes.MN_base_material()

class MN_OT_Add_Custom_Node_Group(bpy.types.Operator):
    bl_idname = "mn.add_custom_node_group"
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
    node_label: bpy.props.StringProperty(name = 'node_label', default = '')
    node_description: bpy.props.StringProperty(
        name = "node_description", 
        description="", 
        default="Add MolecularNodes custom node group.", 
        subtype="NONE"
    )
    node_link: bpy.props.BoolProperty(name = 'node_link', default = True)

    @classmethod
    def poll(cls, context):
        return True
    
    @classmethod
    def description(cls, context, properties):
        return properties.node_description
    
    def execute(self, context):
        try:
            nodes.append(self.node_name, link = self.node_link)
            MN_add_node(self.node_name, label=self.node_label)
        except RuntimeError:
            self.report({'ERROR'}, 
                        message='Failed to add node. Ensure you are not in edit mode.')
        return {"FINISHED"}
    
    def invoke(self, context, event):
        return self.execute(context)

def menu_item_interface(layout_function, 
                        label, 
                        node_name, 
                        node_description='Add custom MolecularNodes node group.', 
                        node_link = True
                        ):
    op = layout_function.operator('mn.add_custom_node_group', text = label)
    op.node_label = label
    op.node_name = node_name
    op.node_description = node_description
    op.node_link = node_link

class MN_OT_Style_Surface_Custom(bpy.types.Operator):
    bl_idname = "mn.style_surface_custom"
    bl_label = "My Class Name"
    bl_description = "Create a split surface representation.\nGenerates an isosurface \
        based on atomic vdw_radii. Each chain has its own separate surface \
        representation"
    bl_options = {"REGISTER", "UNDO"}
    
    @classmethod
    def poll(cls, context):
        return True

    def execute(self, context):
        obj = context.active_object
        try:
            node_surface = nodes.create_custom_surface(
                name = 'MN_style_surface_' + obj.name + '_split', 
                n_chains = len(obj['chain_id_unique'])
            )
        except:
            node_surface = nodes.append('MN_style_surface_single')
            self.report({'WARNING'}, message = 'Unable to detect number of chains.')
        MN_add_node(node_surface.name)
        
        return {"FINISHED"}

class MN_OT_Assembly_Bio(bpy.types.Operator):
    bl_idname = "mn.assembly_bio"
    bl_label = "Build"
    bl_description = "**PDB Downloaded Structures Only**\nAdds node to build \
        biological assembly based on symmetry operations that are extraced from the \
        structure file. Currently this is only supported for structures that were \
        downloaded from the PDB"
    bl_options = {"REGISTER", "UNDO"}
    
    @classmethod
    def poll(cls, context):
        return True

    def execute(self, context):
        from . import assembly
        obj = context.active_object
        
        data_object = assembly.mesh.create_data_object(
            transforms_dict = obj['biological_assemblies'], 
            name = f"data_assembly_{obj.name}"
        )
        
        node_assembly = nodes.create_assembly_node_tree(
            name = obj.name, 
            iter_list = obj['chain_id_unique'], 
            data_object = data_object
            )
        
        MN_add_node(node_assembly.name)
        
        return {"FINISHED"}

def menu_residues_selection_custom(layout_function):
    obj = bpy.context.view_layer.objects.active
    label = 'Res ID'
    op = layout_function.operator(
        'mn.residues_selection_custom', 
        text = label, 
        emboss = True, 
        depress = True
    )

def menu_item_surface_custom(layout_function, label):
    op = layout_function.operator('mn.style_surface_custom', 
                                  text = label, 
                                  emboss = True, 
                                  depress = True)

class MN_OT_Custom_Color_Node(bpy.types.Operator):
    bl_idname = "mn.custom_color_node"
    bl_label = "Custom color by field node."
    bl_options = {"REGISTER", "UNDO"}
    
    node_name: bpy.props.StringProperty(
        name = "node_name", 
        default = ""
    )
    
    node_property: bpy.props.StringProperty(
        name = "node_property", 
        default = ""
    )
    
    field: bpy.props.StringProperty(
        name = "field", 
        default = "chain_id"
    )
    
    prefix: bpy.props.StringProperty(
        name = "prefix", 
        default = "Chain"
    )
    
    @classmethod
    def poll(cls, context):
        return not False
    
    def execute(self, context):
        obj = context.active_object
        try:
            node_color = nodes.chain_color(
                node_name = f"MN_color_{self.node_name}_{obj.name}", 
                input_list = obj[self.node_property], 
                field = self.field, 
                label_prefix= self.prefix
            )
            MN_add_node(node_color.name)
        except:
            self.report({"WARNING"}, message = f"{self.node_propperty} not available for object.")
        return {"FINISHED"}
    
    def invoke(self, context, event):
        return self.execute(context)

def menu_chain_selection_custom(layout_function):
    obj = bpy.context.view_layer.objects.active
    label = 'Chain ' + str(obj.name)
    op = layout_function.operator(
        'mn.chain_selection_custom', 
        text = label, 
        emboss = True, 
        depress = True
    )

class MN_OT_Chain_Selection_Custom(bpy.types.Operator):
    bl_idname = "mn.chain_selection_custom"
    bl_label = "Chain Selection"
    bl_description = "Create a selection based on the chains.\nThis node is built on a \
        per-molecule basis, taking into account the chain_ids that were detected. If \
        no chain information is available this node will not work"
    bl_options = {"REGISTER", "UNDO"}
    
    field: bpy.props.StringProperty(name = "field", default = "chain_id")
    prefix: bpy.props.StringProperty(name = "prefix", default = "Chain ")
    node_property: bpy.props.StringProperty(name = "node_property", default = "chain_id_unique")
    node_name: bpy.props.StringProperty(name = "node_name", default = "chain")
    
    @classmethod
    def poll(cls, context):
        return True
    
    def execute(self, context):
        obj = bpy.context.view_layer.objects.active
        node_chains = nodes.chain_selection(
            node_name = f'MN_select_{self.node_name}_{obj.name}',
            input_list = obj[self.node_property], 
            starting_value = 0,
            attribute = self.field, 
            label_prefix = self.prefix
            )
        
        MN_add_node(node_chains.name)
        
        return {"FINISHED"}


class MN_OT_Residues_Selection_Custom(bpy.types.Operator):
    bl_idname = "mn.residues_selection_custom"
    bl_label = "Multiple Residue Selection"
    bl_description = "Create a selection based on the provided residue strings.\nThis \
        node is built on a per-molecule basis, taking into account the residues that \
        were input."
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
            node_name = 'MN_select_residues', 
            input_resid_string = self.input_resid_string, 
            )
    
        
        MN_add_node(node_residues.name)
        return {"FINISHED"}

    def invoke(self, context, event):
        return context.window_manager.invoke_props_dialog(self)

def menu_ligand_selection_custom(layout_function):
    obj = bpy.context.view_layer.objects.active
    label = 'Ligands ' + str(obj.name)
    op = layout_function.operator(
        'mn.ligand_selection_custom', 
        text = label, 
        emboss = True, 
        depress = True
    )

class MN_OT_Ligand_Selection_Custom(bpy.types.Operator):
    bl_idname = "mn.ligand_selection_custom"
    bl_label = "Ligand Selection"
    bl_description = "Create a selection based on the ligands.\nThis node is built on \
        a per-molecule basis, taking into account the chain_ids that were detected. If \
        no chain information is available this node will not work"
    bl_options = {"REGISTER", "UNDO"}
    
    @classmethod
    def poll(cls, context):
        return True
    
    def execute(self, context):
        obj = bpy.context.view_layer.objects.active
        node_chains = nodes.chain_selection(
            node_name = 'MN_select_' + str(obj.name) + "_ligands", 
            input_list = obj['ligands'], 
            starting_value = 100, 
            attribute = 'res_name', 
            label_prefix = ""
            )
        
        MN_add_node(node_chains.name)
        
        return {"FINISHED"}

class MN_MT_Add_Node_Menu_Properties(bpy.types.Menu):
    bl_idname = 'MN_MT_ADD_NODE_MENU_PROPERTIES'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        # currently nothing for this menu in the panel

class MN_MT_Add_Node_Menu_Color(bpy.types.Menu):
    bl_idname = 'MN_MT_ADD_NODE_MENU_COLOR'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Set Color', 'MN_color_set', 
                            "Sets a new color for the selected atoms")
        layout.separator()
        menu_item_interface(layout, 'Goodsell Colors', 'MN_color_goodsell', 
                            "Adjusts the given colors to copy the 'Goodsell Style'.\n \
                            Darkens the non-carbon atoms and keeps the carbon atoms \
                            the same color. Highlights differences without being too \
                            visually busy")
        layout.separator()
        # menu_item_interface(layout, 'Color by B Factor', 'MN_color_map_attribute')
        menu_item_interface(layout, 'Attribute Map', 'MN_color_attribute_map')
        menu_item_interface(layout, 'Attribute Random', 'MN_color_attribute_random')
        layout.separator()
        menu_item_interface(layout, 'Element', 'MN_color_element', 
                            "Choose a color for each of the first 20 elements")
        # menu_item_color_chains(layout, 'Color by Chains')
        op = layout.operator('mn.custom_color_node', text = 'Chain')
        op.node_property = 'chain_id_unique'
        op.node_name = "chain"
        op.prefix = 'Chain '
        op.field = 'chain_id'
        op = layout.operator('mn.custom_color_node', text = 'Entity')
        op.node_property = 'entity_names'
        op.node_name = "chain"
        op.prefix = ""
        op.field = 'entity_id'
        menu_item_interface(layout, 'Secondary Structure', 'MN_color_sec_struct', 
                            "Specify colors based on the secondary structure")
        menu_item_interface(layout, 'Atomic Number', 'MN_color_atomic_number',
                            "Creates a color based on atomic_number field")
        menu_item_interface(layout, 'Element Common', 'MN_color_common', 
                            "Choose a color for the most common elements in PDB \
                            structures")

class MN_MT_Add_Node_Menu_Bonds(bpy.types.Menu):
    bl_idname = 'MN_MT_ADD_NODE_MENU_BONDS'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Find Bonds', 'MN_bonds_find', 
                            "Finds bonds between atoms based on distance.\n\
                            Based on the vdw_radii for each point, finds other points \
                            within a certain radius to create a bond to. Does not \
                            preserve the index for the points. Does not detect bond type")
        menu_item_interface(layout, 'Break Bonds', 'MN_bonds_break', 
                            "Will delete a bond between atoms that already exists \
                            based on a distance cutoff")
        menu_item_interface(layout, 'Find Bonded Atoms', 'MN_bonds_find_bonded', 
                            "Based on an initial selection, finds atoms which are \
                            within a certain number of bonds away")

class MN_MT_Add_Node_Menu_Styling(bpy.types.Menu):
    bl_idname = 'MN_MT_ADD_NODE_MENU_SYLING'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Atoms', 'MN_style_atoms', 
                            'A sphere atom representation, visible ONLY in Cycles. \
                            Based on point-cloud rendering')
        menu_item_interface(layout, 'Cartoon', 'MN_style_cartoon', 
                            'Create a cartoon representation, highlighting secondary \
                            structure through arrows and ribbons.')
        menu_item_interface(layout, 'Ribbon Protein', 'MN_style_ribbon_protein', 
                            'Create a ribbon mesh based off of the alpha-carbons of \
                            the structure')
        menu_item_interface(layout, 'Ribbon Nucleic', 'MN_style_ribbon_nucleic', 
                            'Create a ribbon mesh and instanced cylinders for nucleic \
                            acids.')
        menu_item_surface_custom(layout, 'Surface Split Chains')
        menu_item_interface(layout, 'Ball and Stick', 'MN_style_ball_and_stick', 
                            "A style node to create ball and stick representation. \
                            Icospheres are instanced on atoms and cylinders for bonds. \
                            Bonds can be detected if they are not present in the \
                            structure")
        menu_item_interface(layout, 'Sticks', 'MN_style_sticks', 
                            "Turn each bond into a cylinder mesh")
        layout.separator()
        layout.label(text = 'Utilities')
        menu_item_interface(layout, 'Atoms Cycles', 'MN_style_atoms_cycles', 
                            'A sphere atom representation, visible ONLY in Cycles. \
                            Based on point-cloud rendering')
        menu_item_interface(layout, 'Atoms EEVEE', 'MN_style_atoms_eevee', 
                            'A sphere atom representation, visible in EEVEE and \
                            Cycles. Based on mesh instancing which slows down viewport \
                            performance')
        menu_item_interface(layout, 'Surface', 'MN_style_surface_single', 
                            "Create a single joined surface representation. \
                            Generates an isosurface based on atomic vdw_radii. All \
                            chains are part of the same surface. Use Surface Split \
                            Chains to have a single surface per chain")
        menu_item_interface(layout, 'Cartoon Utilities', 'MN_style_cartoon_utils')

class MN_MT_Add_Node_Menu_Selections(bpy.types.Menu):
    bl_idname = 'MN_MT_ADD_NODE_MENU_SELECTIONS'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Select Atoms', 'MN_select_atoms', 
                            "Separate atoms based on a selection field.\n" +
                            "Takes atoms and splits them into the selected atoms the \
                            inverted atoms, based on a selection field")
        menu_item_interface(layout, 'Separate Polymers', 'MN_separate_polymers', 
                            "Separate the Geometry into the different polymers.\n" + 
                            "Outputs for protein, nucleic & sugars")
        layout.separator()
        menu_chain_selection_custom(layout)
        op = layout.operator('mn.chain_selection_custom', text = 'Chain')
        op.field = 'chain_id'
        op.prefix = 'Chain '
        op.node_property = 'chain_id_unique'
        op.field = 'chain_id'
        op.node_name = 'chain'
        op = layout.operator('mn.chain_selection_custom', text = 'Entity')
        op.field = 'entity_id'
        op.prefix = ''
        op.node_property = 'entity_names'
        op.field = 'entity_id'
        op.node_name = 'entity'
        menu_ligand_selection_custom(layout)
        menu_item_interface(layout, 'Secondary Structure', 'MN_select_sec_struct')
        layout.separator()
        menu_item_interface(layout, 'Backbone', 'MN_select_backbone', 
                            "Select atoms it they are part of the side chains or backbone.")
        menu_item_interface(layout, 'Atomic Number', 'MN_select_atomic_number', 
                            "Create a selection if input value equal to the \
                            atomic_number field.")
        menu_item_interface(layout, 'Element', 'MN_select_element', 
                            "Create a selection of particular elements by name. Only \
                            first 20 elements supported")
        layout.separator()
        menu_item_interface(layout, 'Proximity', 'MN_select_proximity', 
                            "Select atoms within a certain proximity of some target atoms.")
        menu_item_interface(layout, 'Cube', 'MN_select_cube', 
                            "Create a selection using an Empty Cube", 
                            node_link = False)
        menu_item_interface(layout, 'Sphere', 'MN_select_sphere', 
                            "Create a selection using an Empty Sphere", 
                            node_link = False)
        layout.separator()
        menu_residues_selection_custom(layout)                        
        menu_item_interface(layout, 'Res ID Single', 'MN_select_res_ID_single', 
                            "Create a selection if res_id matches input field")
        menu_item_interface(layout, 'Res ID Range', 'MN_select_res_ID_range', 
                            "Create a selection if the res_id is within the given \
                            thresholds")
        menu_item_interface(layout, 'Res Name Peptide', 'MN_select_res_name', 
                            "Create a selection of particular amino acids by name")
        menu_item_interface(layout, 'Res Name Nucleic', 'MN_select_nucleic_res_name', 
                            "Create a selection of particular nucleic acids by name")
        menu_item_interface(layout, 'Res Whole', 'MN_select_whole_res', 
                            "Expand the selection to every atom in a residue, if any \
                            of those atoms are in the initial selection")

class MN_MT_Add_Node_Menu_Assembly(bpy.types.Menu):
    bl_idname = 'MN_MT_ADD_NODE_MENU_ASSEMBLY'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        layout.operator("mn.assembly_bio", 
                        text = "Biological Assembly", 
                        emboss = True, 
                        depress=True
                        )
        menu_item_interface(layout, 'Center Assembly', 'MN_assembly_center', 
                            "Center the structure on the world origin based on \
                            bounding box")

class MN_MT_Add_Node_Menu_Membranes(bpy.types.Menu):
    bl_idname = 'MN_MT_ADD_NODE_MENU_MEMBRANES'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Setup Atomic Properties', 'MN_prop_setup')

class MN_MT_Add_Node_Menu_DNA(bpy.types.Menu):
    bl_idname = 'MN_MT_ADD_NODE_MENU_DNA'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Double Helix', 'MN_dna_double_helix', 
                            "Create a DNA double helix from an input curve.\n" + 
                            "Takes an input curve and instances for the bases, returns \
                            instances of the bases in a double helix formation")
        menu_item_interface(layout, 'Bases', 'MN_dna_bases', 
                            "Provide the DNA bases as instances to be styled and \
                            passed onto the Double Helix node")
        layout.separator()
        menu_item_interface(layout, 'Style Atoms Cyeles', 'MN_dna_style_atoms_cycles', 
                            "Style the DNA bases with spheres only visible in Cycles")
        menu_item_interface(layout, 'Style Atoms EEVEE', 'MN_dna_style_atoms_eevee', 
                            "Style the DNA bases with spheres visible in Cycles and \
                            EEVEE")
        menu_item_interface(layout, 'Style Surface', 'MN_dna_style_surface', 
                            "Style the DNA bases with surface representation")
        menu_item_interface(layout, 'Style Ball and Stick', 
                            'MN_dna_style_ball_and_stick', 
                            "Style the DNA bases with ball and stick representation")

class MN_MT_Add_Node_Menu_Animation(bpy.types.Menu):
    bl_idname = 'MN_MT_ADD_NODE_MENU_ANIMATION'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Animate Frames', 'MN_animate_frames', 
                            "Interpolate between frames of a trajectory." + 
                            "Given a collection of frames for a trajectory, this node \
                            interpolates between them from start to finish based on \
                            the Animate field taking a value from 0 to 1. The \
                            positions of the Atoms are then moved based on this field")
        menu_item_interface(layout, 'Animate Field', 'MN_animate_field')
        menu_item_interface(layout, 'Animate Value', 'MN_animate_value', 
                            "Animates between given start and end values, based on \
                            the input start and end frame of the timeline. Clamped \
                            will limit the output to the 'To Min' and 'To Max', while \
                            unclamped will continue to interpolate past these values. \
                            'Smoother Step' will ease in and out of these values, with \
                            default being linear interpolation")
        layout.separator()
        menu_item_interface(layout, 'Res Wiggle', "MN_animate_res_wiggle", 
                            "Wiggles the side chains of amino acids based on b_factor, \
                            adding movement to a structure.")
        menu_item_interface(layout, 'Res to Curve', "MN_animate_res_to_curve", 
                            "Takes atoms and maps them along a curve, as a single \
                            long peptide chain.")
        layout.separator()
        menu_item_interface(layout, 'Noise Position', 'MN_noise_position', 
                            "Generate 3D noise field based on the position attribute")
        menu_item_interface(layout, 'Noise Field', 'MN_noise_field', 
                            "Generate a 3D noise field based on the given field")
        menu_item_interface(layout, 'Noise Repeat', 'MN_noise_repeat', 
                            "Generate a 3D noise field that repeats, based on the \
                            given field")

class MN_MT_Add_Node_Menu_Utilities(bpy.types.Menu):
    bl_idname = 'MN_MT_ADD_NODE_MENU_UTILITIES'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Booelean Chain', '.MN_utils_bool_chain')
        menu_item_interface(layout, 'Rotation Matrix', 'MN_utils_rotation_matrix')
        menu_item_interface(layout, 'Curve Resample', 'MN_utils_curve_resample')
        menu_item_interface(layout, 'Determine Secondary Structure', 'MN_utils_dssp')

class MN_MT_Add_Node_Menu_Density(bpy.types.Menu):
    bl_idname = 'MN_MT_ADD_NODE_MENU_DENSITY'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        menu_item_interface(layout, 'Style Surface', 'MN_style_density_surface')
        menu_item_interface(layout, 'Style Wire', 'MN_style_density_wire')
        menu_item_interface(layout, 'Sample Nearest Attribute', 'MN_utils_sample_searest')

class MN_MT_Add_Node_Menu(bpy.types.Menu):
    bl_idname = "MN_MT_ADD_NODE_MENU"
    bl_label = "Menu for Adding Nodes in GN Tree"

    @classmethod
    def poll(cls, context):
        return not (False)

    def draw(self, context):
        layout = self.layout.column_flow(columns=1)
        layout.operator_context = "INVOKE_DEFAULT"
        layout.menu('MN_MT_ADD_NODE_MENU_SYLING', 
                    text='Style', icon_value=77)
        layout.menu('MN_MT_ADD_NODE_MENU_COLOR', 
                    text='Color', icon = 'COLORSET_07_VEC')
        layout.menu('MN_MT_ADD_NODE_MENU_DENSITY', icon = "LIGHTPROBE_CUBEMAP", 
                    text = "Density")
        layout.menu('MN_MT_ADD_NODE_MENU_BONDS', 
                    text='Bonds', icon = 'FIXED_SIZE')
        layout.menu('MN_MT_ADD_NODE_MENU_SELECTIONS', 
                    text='Selection', icon_value=256)
        layout.menu('MN_MT_ADD_NODE_MENU_ANIMATION', 
                    text='Animation', icon_value=409)
        layout.menu('MN_MT_ADD_NODE_MENU_ASSEMBLY', 
                    text='Assemblies', icon = 'GROUP_VERTEX')
        layout.menu('MN_MT_ADD_NODE_MENU_DNA', 
                    text='DNA', icon='GP_SELECT_BETWEEN_STROKES')
        layout.menu('MN_MT_ADD_NODE_MENU_UTILITIES', 
                    text='Utilities', icon_value=92)

def MN_add_node_menu(self, context):
    if ('GeometryNodeTree' == bpy.context.area.spaces[0].tree_type):
        layout = self.layout
        layout.menu('MN_MT_ADD_NODE_MENU', text='Molecular Nodes', icon_value=88)
