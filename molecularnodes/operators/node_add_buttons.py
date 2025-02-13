import bpy
from bpy.types import Context, Operator
from bpy.props import BoolProperty, EnumProperty, IntProperty, StringProperty

from ..blender import nodes
import databpy
from ..ui import node_info


def node_under_mouse(context, event):
    space = context.space_data
    mouse_pos = (event.mouse_region_x, event.mouse_region_y)

    # Find the node under the mouse
    node_under_mouse = None
    for node in space.node_tree.nodes:
        if (
            node.location.x < mouse_pos[0] < node.location.x + node.width
            and node.location.y < mouse_pos[1] < node.location.y + node.height
        ):
            node_under_mouse = node
            break

    return node_under_mouse


def _add_node(node_name, context, show_options=False, material="default"):
    """
    Add a node group to the node tree and set the values.

    intended to be called upon button press in the node tree, and not for use in general scripting
    """

    # actually invoke the operator to add a node to the current node tree
    # use_transform=True ensures it appears where the user's mouse is and is currently
    # being moved so the user can place it where they wish
    bpy.ops.node.add_node(
        "INVOKE_DEFAULT", type="GeometryNodeGroup", use_transform=True
    )
    node = context.active_node
    node.node_tree = bpy.data.node_groups[node_name]
    node.width = nodes.NODE_WIDTH
    node.show_options = show_options
    node.label = node_name
    node.name = node_name

    # if added node has a 'Material' input, set it to the default MN material
    nodes.assign_material(node, new_material=material)


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
            return {"CANCELLED"}
        return {"FINISHED"}


class MN_OT_Assembly_Bio(Operator):
    bl_idname = "mn.assembly_bio"
    bl_label = "Build Biological Assembly"
    bl_description = "Adds node to build biological assembly based on symmetry operations that are extraced from the structure file"
    bl_options = {"REGISTER", "UNDO"}

    inset_node: BoolProperty(default=False)  # type: ignore

    @classmethod
    def poll(self, context):
        # this just checks to see that there is some biological assembly information that
        # is associated with the object / molecule. If there isn't then the assembly
        # operator will be greyed out and unable to be executed
        obj = context.active_object
        if obj is None:
            return False
        return obj.mn.biological_assemblies != ""

    def execute(self, context):
        obj = context.active_object
        with databpy.nodes.DuplicatePrevention():
            try:
                if self.inset_node:
                    nodes.assembly_insert(obj)
                else:
                    tree_assembly = nodes.assembly_initialise(obj)
                    _add_node(tree_assembly.name, context)
            except (KeyError, ValueError) as e:
                self.report({"ERROR"}, "Unable to build biological assembly node.")
                self.report({"ERROR"}, str(e))
                return {"CANCELLED"}

        return {"FINISHED"}


class MN_OT_iswitch_custom(Operator):
    bl_idname = "mn.iswitch_custom"
    # bl_idname = "mn.selection_custom"
    bl_label = "Custom ISwitch Node"
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
    def poll(cls, context: Context) -> bool:
        obj = context.active_object
        if obj is None:
            return False
        return True

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

        prefix = {"BOOLEAN": "Select", "RGBA": "Color"}[self.dtype]
        node_name = " ".join([prefix, self.node_name, name])

        with databpy.nodes.DuplicatePrevention():
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
    bl_label = "Res ID Custom"
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
        with databpy.nodes.DuplicatePrevention():
            node_residues = nodes.resid_multiple_selection(
                node_name="MN_select_res_id_custom",
                input_resid_string=self.input_resid_string,
            )

        _add_node(node_residues.name, context)
        return {"FINISHED"}

    def invoke(self, context, event):
        return context.window_manager.invoke_props_dialog(self)


def get_swap_items(self, context):
    node = context.active_node
    prefix = node.node_tree.name.split(" ")[0].lower()
    if prefix == "is":
        prefix = "select"

    items = [
        (item.name, item.label, item.description)
        for item in node_info.menu_items.get_submenu(prefix).items
        if (not item.is_break and not item.is_custom and item.name != "Set Color")
    ]
    return items


class MN_OT_Node_Swap(Operator):
    bl_idname = "mn.node_swap"
    bl_label = "Swap Node"
    bl_description = "Swap this node for another."

    node_description: StringProperty(default="Swap selected node for another")  # type: ignore
    node_items: EnumProperty(items=get_swap_items)  # type: ignore

    @classmethod
    def description(cls, context, properties):
        return properties.node_description

    def execute(self, context: Context):
        node = context.active_node
        nodes.swap(node, self.node_items)
        return {"FINISHED"}


class MN_OT_Change_Color(Operator):
    bl_idname = "mn.change_color"
    bl_label = "Color"

    color: EnumProperty(  # type: ignore
        items=(
            (item.name, item.label, item.description)
            for item in node_info.menu_items.get_submenu("color").items
            if (not item.is_break and not item.is_custom and item.name != "Set Color")
        )
    )

    def execute(self, context: Context):
        node = context.active_node
        nodes.swap(node, self.color)
        self.report({"INFO"}, f"Selected {self.color}")
        return {"FINISHED"}


CLASSES = [
    MN_OT_Add_Custom_Node_Group,
    MN_OT_Residues_Selection_Custom,
    MN_OT_Assembly_Bio,
    MN_OT_iswitch_custom,
    MN_OT_Change_Color,
    MN_OT_Node_Swap,
]
