import bpy
from ..blender import nodes





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
            nodes.add_node(self.node_name, label=self.node_label)
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
                        node_link = False
                        ):
    op = layout_function.operator('mn.add_custom_node_group', text = label)
    op.node_label = label
    op.node_name = node_name
    op.node_description = node_description
    op.node_link = node_link

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
        transforms_array = assembly.mesh.get_transforms_from_dict(obj['biological_assemblies'])
        data_object = assembly.mesh.create_data_object(
            transforms_array = transforms_array, 
            name = f"data_assembly_{obj.name}"
        )
        
        node_assembly = nodes.create_assembly_node_tree(
            name = obj.name, 
            iter_list = obj['chain_id_unique'], 
            data_object = data_object
            )
        
        nodes.add_node(node_assembly.name)
        
        return {"FINISHED"}

class MN_OT_Color_Custom(bpy.types.Operator):
    bl_idname = "mn.color_custom"
    bl_label = "Custom color by field node."
    bl_options = {"REGISTER", "UNDO"}
    
    description: bpy.props.StringProperty(name = "description", default = "")
    
    node_name: bpy.props.StringProperty(name = "node_name", default = "")
    node_property: bpy.props.StringProperty(name = "node_property", default = "")
    field: bpy.props.StringProperty(name = "field", default = "chain_id")
    prefix: bpy.props.StringProperty(name = "prefix", default = "Chain")
    starting_value: bpy.props.IntProperty(name = "starting_value", default = 0)
    
    @classmethod
    def poll(cls, context):
        return not False
    
    @classmethod
    def description(cls, context, properties):
        return properties.descriptionf
    
    def execute(self, context):
        obj = context.active_object
        # try:
        node_color = nodes.chain_color(
            node_name = f"MN_color_{self.node_name}_{obj.name}", 
            input_list = obj[self.node_property], 
            field = self.field, 
            label_prefix= self.prefix, 
            starting_value = self.starting_value
        )
        nodes.add_node(node_color.name)
        # except:
            # self.report({"WARNING"}, message = f"{self.node_propperty} not available for object.")
        return {"FINISHED"}
    
    def invoke(self, context, event):
        return self.execute(context)

class MN_OT_selection_custom(bpy.types.Operator):
    bl_idname = "mn.selection_custom"
    bl_label = "Chain Selection"
    bl_options = {"REGISTER", "UNDO"}
    
    description: bpy.props.StringProperty(name = "description", default = "")
    
    field: bpy.props.StringProperty(name = "field", default = "chain_id")
    prefix: bpy.props.StringProperty(name = "prefix", default = "Chain ")
    node_property: bpy.props.StringProperty(name = "node_property", default = "chain_id_unique")
    node_name: bpy.props.StringProperty(name = "node_name", default = "chain")
    starting_value: bpy.props.IntProperty(name = "starting_value", default = 0)
    
    @classmethod
    def poll(cls, context):
        return True
    
    @classmethod
    def description(cls, context, properties):
        return properties.description
    
    def execute(self, context):
        obj = context.view_layer.objects.active
        node_chains = nodes.chain_selection(
            node_name = f'MN_select_{self.node_name}_{obj.name}',
            input_list = obj[self.node_property], 
            starting_value = self.starting_value,
            attribute = self.field, 
            label_prefix = self.prefix
            )
        
        nodes.add_node(node_chains.name)
        
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
            node_name = 'MN_select_res_id_custom', 
            input_resid_string = self.input_resid_string, 
            )
    
        nodes.add_node(node_residues.name)
        return {"FINISHED"}

    def invoke(self, context, event):
        return context.window_manager.invoke_props_dialog(self)

def button_custom_color(layout, text, field, prefix, property, node_name, starting_value = 0):
    op = layout.operator('mn.color_custom', text = text)
    op.field = field
    op.prefix = prefix
    op.node_property = property
    op.node_name = node_name
    op.starting_value = starting_value
    op.description = f"Choose individual colors for each {text}"

def button_custom_selection(layout, text, field, prefix, property, node_name, starting_value = 0):
    op = layout.operator('mn.selection_custom', text = text)
    op.field = field
    op.prefix = prefix
    op.node_property = property
    op.node_name = node_name
    op.starting_value = starting_value
    op.description = f"Create individual selections for each {text}"

