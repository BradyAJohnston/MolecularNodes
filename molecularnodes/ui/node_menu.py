import bpy

from .node_info import menu_items
from .func import build_menu


class MN_MT_Node_Color(bpy.types.Menu):
    bl_idname = 'MN_MT_NODE_COLOR'
    bl_label = ''
    
    def draw(self, context):
        layout = self.layout
        build_menu(layout, menu_items['color'])


class MN_MT_Node_Bonds(bpy.types.Menu):
    bl_idname = 'MN_MT_NODE_BONDS'
    bl_label = ''
    
    def draw(self, context):
        layout = self.layout
        build_menu(layout, menu_items['bonds'])


class MN_MT_Node_Style(bpy.types.Menu):
    bl_idname = 'MN_MT_NODE_STYLE'
    bl_label = ''
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        build_menu(layout, menu_items['style'])

class MN_MT_Node_Select(bpy.types.Menu):
    bl_idname = 'MN_MT_NODE_SELECT'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        build_menu(layout, menu_items['select'])

class MN_MT_Node_Assembly(bpy.types.Menu):
    bl_idname = 'MN_MT_NODE_ASSEMBLY'
    bl_label = ''
    
    @classmethod
    def poll(cls, context):
        return True
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        build_menu(layout, menu_items['assembly'])

class MN_MT_Node_DNA(bpy.types.Menu):
    bl_idname = 'MN_MT_NODE_DNA'
    bl_label = ''
    
    def draw(self, context):
        layout = self.layout
        build_menu(layout, menu_items['dna'])

class MN_MT_Node_Animate(bpy.types.Menu):
    bl_idname = 'MN_MT_NODE_ANIMATE'
    bl_label = ''
    
    def draw(self, context):
        layout = self.layout
        build_menu(layout, menu_items['animate'])

class MN_MT_Node_Utils(bpy.types.Menu):
    bl_idname = 'MN_MT_NODE_UTILS'
    bl_label = ''
    
    def draw(self, context):
        build_menu(self.layout, menu_items['utils'])

class MN_MT_Node_CellPack(bpy.types.Menu):
    bl_idname = "MN_MT_NODE_CELLPACK"
    bl_label = ''
    
    def draw(self, context):
        layout = self.layout
        build_menu(layout, menu_items['cellpack'])

class MN_MT_Node_Density(bpy.types.Menu):
    bl_idname = 'MN_MT_NODE_DENSITY'
    bl_label = ''
    
    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_DEFAULT"
        build_menu(layout, menu_items['density'])


class MN_MT_Node(bpy.types.Menu):
    bl_idname = "MN_MT_NODE"
    bl_label = "Menu for Adding Nodes in GN Tree"

    @classmethod
    def poll(cls, context):
        return True

    def draw(self, context):
        layout = self.layout.column_flow(columns=1)
        layout.operator_context = "INVOKE_DEFAULT"

        layout.menu('MN_MT_NODE_STYLE', text='Style', icon_value=77)
        layout.menu('MN_MT_NODE_SELECT', text='Selection', icon_value=256)
        layout.menu('MN_MT_NODE_COLOR', text='Color', icon='COLORSET_07_VEC')
        layout.menu('MN_MT_NODE_ANIMATE', text='Animation', icon_value=409)
        layout.menu('MN_MT_NODE_ASSEMBLY', text='Assemblies', icon='GROUP_VERTEX')
        layout.menu('MN_MT_NODE_CELLPACK', text='CellPack', icon='PARTICLE_POINT')
        layout.menu('MN_MT_NODE_DENSITY', text='Density', icon="VOLUME_DATA")
        layout.menu('MN_MT_NODE_DNA', text='DNA', icon='GP_SELECT_BETWEEN_STROKES')
        layout.menu('MN_MT_NODE_UTILS', text='Utilities', icon_value=92)

def MN_add_node_menu(self, context):
    if ('GeometryNodeTree' == bpy.context.area.spaces[0].tree_type):
        layout = self.layout
        layout.menu('MN_MT_NODE', text='Molecular Nodes', icon_value=88)
