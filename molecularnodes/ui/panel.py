import bpy
from .. import pkg
from ..io import (
    pdb, local, star, cellpack, md, density
)

bpy.types.Scene.MN_panel_import = bpy.props.EnumProperty(
    name = "Import Method", 
    items = (
        ('pdb', "PDB", "Download from the PDB", 0),
        ('local', "Local", "Open a local file", 1),
        ('md', "MD", "Import a molecular dynamics trajectory", 2),
        ('density', "Density", "Import an EM Density Map", 3), 
        ('star', 'Starfile', "Import a .starfile mapback file", 4), 
        ('cellpack', 'CellPack', "Import a CellPack .cif/.bcif file.", 5)
    )
)
bpy.types.Scene.MN_panel = bpy.props.EnumProperty(
    name = "Import Method", 
    items = (
        ('import', "Import", "Import macromolecules", 0),
        ('object', "Object", "Adjust settings affecting the selected object.", 1),
        ('scene',  "Scene", "Change settings for the world and rendering", 2)
    )
)

chosen_panel = {
    'pdb': pdb, 
    'local': local, 
    'star': star, 
    'md': md, 
    'density': density, 
    'cellpack': cellpack, 
    
}

packages = {
    'pdb': ['biotite', 'scipy'], 
    'star': ['starfile', 'eulerangles'], 
    'local': ['biotite', 'scipy'], 
    'cellpack': ['biotite', 'msgpack', 'scipy'], 
    'md': ['MDAnalysis'], 
    'density': ['mrcfile', 'scipy'], 
}


def panel_import(layout, scene):
    layout.label(text = "Import Options", icon = "MODIFIER")
    box = layout.box()
    grid = box.grid_flow(columns = 2)
    
    grid.prop(scene, 'MN_import_centre', icon_value=0)
    grid.prop(scene, 'MN_import_del_solvent', icon_value=0)
    grid.prop(scene, "MN_import_style")
    row = layout.grid_flow(row_major = True, columns = 3, align = True)
    row.alignment = 'EXPAND'
    row.enabled = True
    row.alert = False
    
    row.prop(scene, 'MN_panel_import')
    
    selection = scene.MN_panel_import
    col = layout.column()
    box = col.box()
    
    row = layout.row()
    for package in packages[selection]:
        if not pkg.is_current(package):
            box.enabled = False
            box.alert = True
            box.label(text = f'Please install {package} in the Molecular Nodes preferences.')
    chosen_panel[selection].panel(box, scene)

def panel_scene(layout, scene):
    layout.label(text = 'adjust scene properties')
    layout.prop(scene.render, 'film_transparent')
    layout.prop(bpy.data.worlds["World Shader"].node_tree.nodes["MN_world_shader"].inputs[2], 'default_value')
# node_tree.nodes["MN_world_shader"].inputs[2].default_value
class MN_PT_panel(bpy.types.Panel):
    bl_label = 'Molecular Nodes'
    bl_idname = 'MN_PT_panel'
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = 'scene'
    bl_order = 0
    bl_options = {'HEADER_LAYOUT_EXPAND'}
    bl_ui_units_x=0

    def draw(self, context):
        layout = self.layout
        scene = context.scene
        layout.prop_tabs_enum(scene, 'MN_panel')
        if scene.MN_panel == "import":
            panel_import(layout, scene)
        elif scene.MN_panel == "scene":
            panel_scene(layout, scene)