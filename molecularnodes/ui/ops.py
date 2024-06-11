import bpy
from ..blender import nodes


class MN_OT_Add_Custom_Node_Group(bpy.types.Operator):
    bl_idname = "mn.add_custom_node_group"
    bl_label = "Add Custom Node Group"
    # bl_description = "Add Molecular Nodes custom node group."
    bl_options = {"REGISTER", "UNDO"}
    node_name: bpy.props.StringProperty(
        name="node_name", description="", default="", subtype="NONE", maxlen=0
    )
    node_label: bpy.props.StringProperty(name="node_label", default="")
    node_description: bpy.props.StringProperty(
        name="node_description",
        description="",
        default="Add MolecularNodes custom node group.",
        subtype="NONE",
    )
    node_link: bpy.props.BoolProperty(name="node_link", default=True)

    @classmethod
    def description(cls, context, properties):
        return properties.node_description

    def execute(self, context):
        try:
            nodes.append(self.node_name, link=self.node_link)
            nodes.add_node(self.node_name)  # , label=self.node_label)
        except RuntimeError:
            self.report(
                {"ERROR"},
                message="Failed to add node. Ensure you are not in edit mode.",
            )
        return {"FINISHED"}


class MN_OT_Assembly_Bio(bpy.types.Operator):
    bl_idname = "mn.assembly_bio"
    bl_label = "Build"
    bl_description = "Adds node to build biological assembly based on symmetry operations that are extraced from the structure file"
    bl_options = {"REGISTER", "UNDO"}

    inset_node: bpy.props.BoolProperty(default=False)

    @classmethod
    def poll(self, context):
        # this just checks to see that there is some biological assembly information that
        # is associated with the object / molecule. If there isn't then the assembly
        # operator will be greyed out and unable to be executed
        bob = context.active_object
        try:
            bob["biological_assemblies"]
            return True
        except KeyError:
            False

    def execute(self, context):
        bob = context.active_object
        try:
            if self.inset_node:
                nodes.assembly_insert(bob)
            else:
                tree_assembly = nodes.assembly_initialise(bob)
                nodes.add_node(tree_assembly.name)
        except (KeyError, ValueError) as e:
            self.report({"ERROR"}, "Unable to build biological assembly node.")
            self.report({"ERROR"}, str(e))
            return {"CANCELLED"}

        return {"FINISHED"}


class MN_OT_iswitch_custom(bpy.types.Operator):
    bl_idname = "mn.iswitch_custom"
    # bl_idname = "mn.selection_custom"
    bl_label = "Chain Selection"
    bl_options = {"REGISTER", "UNDO"}

    description: bpy.props.StringProperty(name="Description")
    dtype: bpy.props.EnumProperty(  # type: ignore
        name="Data type",
        items=(
            ("RGBA", "RGBA", "Color iswitch."),
            ("BOOLEAN", "BOOLEAN", "Boolean iswitch"),
        ),
    )
    field: bpy.props.StringProperty(name="field", default="chain_id")
    prefix: bpy.props.StringProperty(name="prefix", default="Chain ")
    node_property: bpy.props.StringProperty(name="node_property", default="chain_ids")
    node_name: bpy.props.StringProperty(name="node_name", default="chain")
    starting_value: bpy.props.IntProperty(name="starting_value", default=0)

    @classmethod
    def description(cls, context, properties):
        return properties.description

    def execute(self, context):
        object = context.view_layer.objects.active
        prop = object[self.node_property]
        name = object.name
        if not prop:
            self.report(
                {"WARNING"},
                message=f"{self.node_property} not available for {object.name}.",
            )
            return {"CANCELLED"}

        if self.dtype == "BOOLEAN":
            node_name = f"Select {self.node_name}_{name}"
        elif self.dtype == "RGBA":
            node_name = f"Select {self.node_name}_{name}"
        else:
            raise ValueError(f"Data type not supported {self.dtype}")

        node_chains = nodes.custom_iswitch(
            name=node_name,
            dtype=self.dtype,
            iter_list=prop,
            start=self.starting_value,
            field=self.field,
            prefix=self.prefix,
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
        default="19,94,1-16",
    )

    def execute(self, context):
        node_residues = nodes.resid_multiple_selection(
            node_name="MN_select_res_id_custom",
            input_resid_string=self.input_resid_string,
        )

        nodes.add_node(node_residues.name)
        return {"FINISHED"}

    def invoke(self, context, event):
        return context.window_manager.invoke_props_dialog(self)


class MN_OT_Change_Style(bpy.types.Operator):
    bl_idname = "mn.style_change"
    bl_label = "Style"

    style: bpy.props.EnumProperty(name="Style", items=nodes.STYLE_ITEMS)

    def execute(self, context):
        object = context.active_object
        nodes.change_style_node(object, self.style)
        return {"FINISHED"}


ops_ui = [
    MN_OT_Add_Custom_Node_Group,
    MN_OT_Residues_Selection_Custom,
    MN_OT_Change_Style,
]
