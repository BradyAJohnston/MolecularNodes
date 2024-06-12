import bpy
from bpy.types import Operator
from bpy.props import BoolProperty, EnumProperty, IntProperty, StringProperty

from ..blender import nodes
from .panel import STYLE_ITEMS


def _add_node(
    node_name, context, label: str = "", show_options=False, material="default"
):
    # intended to be called upon button press in the node tree, and not for use
    # in general scripting

    prev_context = context.area.type
    context.area.type = "NODE_EDITOR"
    # actually invoke the operator to add a node to the current node tree
    # use_transform=True ensures it appears where the user's mouse is and is currently
    # being moved so the user can place it where they wish
    bpy.ops.node.add_node(
        "INVOKE_DEFAULT", type="GeometryNodeGroup", use_transform=True
    )
    context.area.type = prev_context
    node = context.active_node
    node.node_tree = bpy.data.node_groups[node_name]
    node.width = 200.0
    node.show_options = show_options

    # if label == "":
    #     node.label = format_node_name(node_name)
    # else:
    #     node.label = label
    node.label = node_name
    node.name = node_name

    # if added node has a 'Material' input, set it to the default MN material
    nodes.assign_material(node, material=material)


class MN_OT_Add_Custom_Node_Group(Operator):
    bl_idname = "mn.add_custom_node_group"
    bl_label = "Add Custom Node Group"
    # bl_description = "Add Molecular Nodes custom node group."
    bl_options = {"REGISTER", "UNDO"}
    node_name: StringProperty(  # type: ignore
        name="node_name", description="", default="", subtype="NONE", maxlen=0
    )
    node_label: StringProperty(name="node_label", default="")  # type: ignore
    node_description: StringProperty(  # type: ignore
        name="node_description",
        description="",
        default="Add MolecularNodes custom node group.",
        subtype="NONE",
    )
    node_link: BoolProperty(name="node_link", default=True)  # type: ignore

    @classmethod
    def description(cls, context, properties):
        return properties.node_description

    def execute(self, context):
        try:
            nodes.append(self.node_name, link=self.node_link)
            _add_node(self.node_name, context)  # , label=self.node_label)
        except RuntimeError:
            self.report(
                {"ERROR"},
                message="Failed to add node. Ensure you are not in edit mode.",
            )
        return {"FINISHED"}


class MN_OT_Assembly_Bio(Operator):
    bl_idname = "mn.assembly_bio"
    bl_label = "Build"
    bl_description = "Adds node to build biological assembly based on symmetry operations that are extraced from the structure file"
    bl_options = {"REGISTER", "UNDO"}

    inset_node: BoolProperty(default=False)  # type: ignore

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
                _add_node(tree_assembly.name, context)
        except (KeyError, ValueError) as e:
            self.report({"ERROR"}, "Unable to build biological assembly node.")
            self.report({"ERROR"}, str(e))
            return {"CANCELLED"}

        return {"FINISHED"}


class MN_OT_iswitch_custom(Operator):
    bl_idname = "mn.iswitch_custom"
    # bl_idname = "mn.selection_custom"
    bl_label = "Chain Selection"
    bl_options = {"REGISTER", "UNDO"}

    description: StringProperty(name="Description")  # type: ignore
    dtype: EnumProperty(  # type: ignore
        name="Data type",
        items=(
            ("RGBA", "RGBA", "Color iswitch."),
            ("BOOLEAN", "BOOLEAN", "Boolean iswitch"),
        ),
    )
    field: StringProperty(name="field", default="chain_id")  # type: ignore
    prefix: StringProperty(name="prefix", default="Chain ")  # type: ignore
    node_property: StringProperty(name="node_property", default="chain_ids")  # type: ignore
    node_name: StringProperty(name="node_name", default="chain")  # type: ignore
    starting_value: IntProperty(name="starting_value", default=0)  # type: ignore

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

        _add_node(node_chains.name, context)

        return {"FINISHED"}


class MN_OT_Residues_Selection_Custom(Operator):
    bl_idname = "mn.residues_selection_custom"
    bl_label = "Multiple Residue Selection"
    bl_description = "Create a selection based on the provided residue strings.\nThis \
        node is built on a per-molecule basis, taking into account the residues that \
        were input."
    bl_options = {"REGISTER", "UNDO"}

    input_resid_string: StringProperty(  # type: ignore
        name="Select residue IDs: ",
        description="Enter a string value.",
        default="19,94,1-16",
    )

    def execute(self, context):
        node_residues = nodes.resid_multiple_selection(
            node_name="MN_select_res_id_custom",
            input_resid_string=self.input_resid_string,
        )

        _add_node(node_residues.name, context)
        return {"FINISHED"}

    def invoke(self, context, event):
        return context.window_manager.invoke_props_dialog(self)


class MN_OT_Change_Style(Operator):
    bl_idname = "mn.style_change"
    bl_label = "Style"

    style: EnumProperty(name="Style", items=STYLE_ITEMS)  # type: ignore

    def execute(self, context):
        object = context.active_object
        nodes.change_style_node(object, self.style)
        return {"FINISHED"}


class MN_OT_Swap_Style_Node(bpy.types.Operator):
    bl_idname = "mn.style_change_node"
    bl_label = "Style"

    style: bpy.props.EnumProperty(name="Style", items=STYLE_ITEMS)  # type: ignore

    @classmethod
    def poll(self, context):
        node = context.space_data.edit_tree.nodes.active
        return node.name.startswith("Style")

    def execute(self, context):
        nodes.swap_style_node(
            tree=context.space_data.node_tree,
            node_style=context.space_data.edit_tree.nodes.active,
            style=self.style,
        )
        return {"FINISHED"}


CLASSES = [
    MN_OT_Add_Custom_Node_Group,
    MN_OT_Residues_Selection_Custom,
    MN_OT_Change_Style,
    MN_OT_Assembly_Bio,
    MN_OT_iswitch_custom,
    MN_OT_Swap_Style_Node,
]
