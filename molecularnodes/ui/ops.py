import bpy
from .. import assembly
from ..blender import nodes

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
    def description(cls, context, properties):
        return properties.node_description
    
    def execute(self, context):
        try:
            nodes.append(self.node_name, link = self.node_link)
            nodes.add_node(self.node_name) #, label=self.node_label)
        except RuntimeError:
            self.report({'ERROR'}, 
                        message='Failed to add node. Ensure you are not in edit mode.')
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
    def poll(self, context):
        mol = context.active_object
        return mol.mn['molecule_type'] in ['pdb', 'local']

    def execute(self, context):
        tree_assembly = nodes.assembly_initialise(context.active_object)
        nodes.add_node(tree_assembly.name)
        
        return {"FINISHED"}


class MN_OT_Color_Custom(bpy.types.Operator):
    bl_idname = "mn.color_custom"
    bl_label = "Custom color by field node."
    bl_options = {"REGISTER", "UNDO"}
    
    description: bpy.props.StringProperty(name = "description", default = "")
    
    node_name: bpy.props.StringProperty(name = "node_name", default = "")
    node_property: bpy.props.StringProperty(name = "node_property", default = "chain_id_unique")
    field: bpy.props.StringProperty(name = "field", default = "chain_id")
    prefix: bpy.props.StringProperty(name = "prefix", default = "Chain")
    starting_value: bpy.props.IntProperty(name = "starting_value", default = 0)
    
    @classmethod
    def description(cls, context, properties):
        return properties.description
    
    def execute(self, context):
        obj = context.active_object
        # try:
        node_color = nodes.chain_color(
            name = f"MN_color_{self.node_name}_{obj.name}", 
            input_list = obj[self.node_property], 
            field = self.field, 
            label_prefix= self.prefix, 
            starting_value = self.starting_value
        )
        nodes.add_node(node_color.name)
        # except:
            # self.report({"WARNING"}, message = f"{self.node_propperty} not available for object.")
        return {"FINISHED"}


class MN_OT_selection_custom(bpy.types.Operator):
    bl_idname = "mn.selection_custom"
    bl_label = "Chain Selection"
    bl_options = {"REGISTER", "UNDO"}
    
    
    description: bpy.props.StringProperty(name = "Description")
    field: bpy.props.StringProperty(name = "field", default = "chain_id")
    prefix: bpy.props.StringProperty(name = "prefix", default = "Chain ")
    node_property: bpy.props.StringProperty(name = "node_property", default = "chain_id_unique")
    node_name: bpy.props.StringProperty(name = "node_name", default = "chain")
    starting_value: bpy.props.IntProperty(name = "starting_value", default = 0)
    
    @classmethod
    def description(cls, context, properties):
        return properties.description
    
    def execute(self, context):
        obj = context.view_layer.objects.active
        node_chains = nodes.chain_selection(
            name = f'MN_select_{self.node_name}_{obj.name}',
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
    
    def execute(self, context):
        node_residues = nodes.resid_multiple_selection(
            node_name = 'MN_select_res_id_custom', 
            input_resid_string = self.input_resid_string, 
            )
    
        nodes.add_node(node_residues.name)
        return {"FINISHED"}
    
    def invoke(self, context, event):
        return context.window_manager.invoke_props_dialog(self)
