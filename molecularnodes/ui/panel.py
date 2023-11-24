import bpy
from .. import pkg
from ..io import (
    pdb, local, star, cellpack, md, density, dna
)

bpy.types.Scene.MN_panel = bpy.props.EnumProperty(
    name = "Panel Selection", 
    items = (
        ('import', "Import", "Import macromolecules", 0),
        ('object', "Object", "Adjust settings affecting the selected object.", 1),
        ('scene',  "Scene", "Change settings for the world and rendering", 2)
    )
)
bpy.types.Scene.MN_panel_import = bpy.props.EnumProperty(
    name = "Import Method", 
    items = (
        ('pdb', "PDB", "Download from the PDB", 0),
        ('local', "Local", "Open a local file", 1),
        ('md', "MD", "Import a molecular dynamics trajectory", 2),
        ('density', "Density", "Import an EM Density Map", 3), 
        ('star', 'Starfile', "Import a .starfile mapback file", 4), 
        ('cellpack', 'CellPack', "Import a CellPack .cif/.bcif file.", 5), 
        ('dna', 'oxDNA', 'Import an oxDNA fil.', 6)
    )
)

chosen_panel = {
    'pdb': pdb, 
    'local': local, 
    'star': star, 
    'md': md, 
    'density': density, 
    'cellpack': cellpack,
    'dna': dna
    
}

packages = {
    'pdb': ['biotite', 'scipy'], 
    'star': ['starfile', 'eulerangles'], 
    'local': ['biotite', 'scipy'], 
    'cellpack': ['biotite', 'msgpack', 'scipy'], 
    'md': ['MDAnalysis'], 
    'density': ['mrcfile', 'scipy'],
    'dna': []
}


def panel_import(layout, scene):
    selection = scene.MN_panel_import
    layout.prop(scene, 'MN_panel_import')
    buttons = layout.grid_flow()
    col = layout.column()
    for package in packages[scene.MN_panel_import]:
        if not pkg.is_current(package):
            pkg.button_install_pkg(buttons, package, pkg.get_pkgs()[package]['version'])
    box = col.box()
    for package in packages[selection]:
        if not pkg.is_current(package):
            box.enabled = False
            box.alert = True
            box.label(text = f'Please install {package} in the Molecular Nodes preferences.')
    chosen_panel[selection].panel(box, scene)


def panel_scene(layout, scene):
    
    cam = bpy.data.cameras[bpy.data.scenes["Scene"].camera.name]
    world_shader = bpy.data.worlds["World Shader"].node_tree.nodes["MN_world_shader"]
    grid = layout.grid_flow()
    world = grid.box().column(heading="World Settings")
    world.prop(bpy.data.scenes["Scene"].render, "engine")
    world.prop(world_shader.inputs[1], 'default_value', text = "World Lighting")
    world.label(text = "Background")
    row = world.row()
    row.prop(scene.render, 'film_transparent')
    row.prop(world_shader.inputs[2], 'default_value', text = "")
    
    camera = grid.box().column(heading="Camera Settings")
    camera.prop(cam, "lens")
    camera.prop(cam.dof, 'use_dof')
    focus = camera.column()
    focus.enabled = cam.dof.use_dof
    focus.prop(cam.dof, 'focus_object')
    distance = focus.row()
    distance.enabled = (cam.dof.focus_object is None)
    distance.prop(cam.dof, 'focus_distance')
    focus.prop(cam.dof, 'aperture_fstop')
    camera.prop(bpy.data.scenes["Scene"].render, "use_motion_blur")


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