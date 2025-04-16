import bpy
from .node_info import menu_items


class MN_MT_Node_Color(bpy.types.Menu):
    bl_idname = "MN_MT_NODE_COLOR"
    bl_label = ""

    def draw(self, context):
        menu_items.get_submenu("color").menu(self.layout, context)


class MN_MT_Node_Bonds(bpy.types.Menu):
    bl_idname = "MN_MT_NODE_BONDS"
    bl_label = ""

    def draw(self, context):
        menu_items.get_submenu("bonds").menu(self.layout, context)


class MN_MT_Node_Style(bpy.types.Menu):
    bl_idname = "MN_MT_NODE_STYLE"
    bl_label = ""

    def draw(self, context):
        menu_items.get_submenu("style").menu(self.layout, context)


class MN_MT_Node_Select(bpy.types.Menu):
    bl_idname = "MN_MT_NODE_SELECT"
    bl_label = ""

    def draw(self, context):
        menu_items.get_submenu("select").menu(self.layout, context)


class MN_MT_Node_Assembly(bpy.types.Menu):
    bl_idname = "MN_MT_NODE_ASSEMBLY"
    bl_label = ""

    def draw(self, context):
        menu_items.get_submenu("ensemble").menu(self.layout, context)


class MN_MT_Node_DNA(bpy.types.Menu):
    bl_idname = "MN_MT_NODE_DNA"
    bl_label = ""

    def draw(self, context):
        menu_items.get_submenu("DNA").menu(self.layout, context)


class MN_MT_Node_Animate(bpy.types.Menu):
    bl_idname = "MN_MT_NODE_ANIMATE"
    bl_label = ""

    def draw(self, context):
        menu_items.get_submenu("animate").menu(self.layout, context)


class MN_MT_Node_Density(bpy.types.Menu):
    bl_idname = "MN_MT_NODE_DENSITY"
    bl_label = ""

    def draw(self, context):
        menu_items.get_submenu("density").menu(self.layout, context)


class MN_MT_Node_Topology(bpy.types.Menu):
    bl_idname = "MN_MT_NODE_TOPOLOGY"
    bl_label = ""

    def draw(self, context):
        menu_items.get_submenu("topology").menu(self.layout, context)


class MN_MT_Node_Attributes(bpy.types.Menu):
    bl_idname = "MN_MT_NODE_ATTRIBUTES"
    bl_label = ""

    def draw(self, context):
        menu_items.get_submenu("attributes").menu(self.layout, context)


class MN_MT_Node_Curves(bpy.types.Menu):
    bl_idname = "MN_MT_NODE_CURVES"
    bl_label = ""

    def draw(self, context):
        menu_items.get_submenu("curves").menu(self.layout, context)


class MN_MT_Node_Geometry(bpy.types.Menu):
    bl_idname = "MN_MT_NODE_GEOMETRY"
    bl_label = ""

    def draw(self, context):
        menu_items.get_submenu("geometry").menu(self.layout, context)


class MN_MT_Node_Simulation(bpy.types.Menu):
    bl_idname = "MN_MT_NODE_SIMULATION"
    bl_label = ""

    def draw(self, context):
        menu_items.get_submenu("simulation").menu(self.layout, context)


class MN_MT_Node_Utils(bpy.types.Menu):
    bl_idname = "MN_MT_NODE_UTILS"
    bl_label = ""

    def draw(self, context):
        menu_items.get_submenu("utils").menu(self.layout, context)


class MN_MT_Node_Fields(bpy.types.Menu):
    bl_idname = "MN_MT_NODE_FIELDS"
    bl_label = ""

    def draw(self, context):
        menu_items.get_submenu("fields").menu(self.layout, context)


def draw_node_menus(self, context):
    layout = self.layout
    layout.menu("MN_MT_NODE_STYLE", text="Style", icon_value=77)
    layout.menu("MN_MT_NODE_SELECT", text="Select", icon="RESTRICT_SELECT_OFF")
    layout.menu("MN_MT_NODE_COLOR", text="Color", icon="COLORSET_07_VEC")
    layout.separator()
    layout.menu("MN_MT_NODE_ANIMATE", text="Animation", icon="MOD_DASH")
    layout.menu("MN_MT_NODE_GEOMETRY", text="Geometry", icon="MESH_DATA")
    layout.menu("MN_MT_NODE_SIMULATION", text="Simulation", icon="PHYSICS")
    layout.menu("MN_MT_NODE_ASSEMBLY", text="Ensemble", icon="GROUP_VERTEX")
    layout.menu("MN_MT_NODE_TOPOLOGY", text="Topology", icon="ORIENTATION_CURSOR")
    layout.menu("MN_MT_NODE_ATTRIBUTES", text="Attributes", icon="SPREADSHEET")
    layout.separator()
    layout.menu("MN_MT_NODE_DENSITY", text="Density", icon="VOLUME_DATA")
    layout.separator()
    layout.menu("MN_MT_NODE_DNA", text="DNA", icon="GP_SELECT_BETWEEN_STROKES")
    layout.separator()
    layout.menu("MN_MT_NODE_CURVES", text="Curves", icon="CURVE_DATA")
    layout.menu("MN_MT_NODE_UTILS", text="Utilities", icon="TOOL_SETTINGS")
    layout.menu("MN_MT_NODE_FIELDS", text="Fields", icon="NODE")


class MN_MT_Node(bpy.types.Menu):
    bl_idname = "MN_MT_NODE"
    bl_label = "Menu for Adding Nodes in GN Tree"

    def draw(self, context):
        draw_node_menus(self, context)


def add_node_menu(self, context):
    if "GeometryNodeTree" == bpy.context.area.spaces[0].tree_type:
        layout = self.layout
        layout.menu("MN_MT_NODE", text="Molecular Nodes", icon="PARTICLE_DATA")


CLASSES = [
    MN_MT_Node,
    MN_MT_Node_Animate,
    MN_MT_Node_Assembly,
    MN_MT_Node_Bonds,
    MN_MT_Node_Color,
    MN_MT_Node_Density,
    MN_MT_Node_DNA,
    MN_MT_Node_Style,
    MN_MT_Node_Select,
    MN_MT_Node_Topology,
    MN_MT_Node_Attributes,
    MN_MT_Node_Geometry,
    MN_MT_Node_Simulation,
    MN_MT_Node_Curves,
    MN_MT_Node_Utils,
    MN_MT_Node_Fields,
]
