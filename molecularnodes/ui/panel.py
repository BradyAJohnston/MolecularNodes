import bpy
from .. import pkg
from ..ui import pref
from ..blender import nodes
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
    name = "Method", 
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
    'pdb': ['biotite'], 
    'star': ['starfile'], 
    'local': ['biotite'], 
    'cellpack': ['biotite', 'msgpack'], 
    'md': ['MDAnalysis'], 
    'density': ['mrcfile'],
    'dna': []
}

class MN_OT_Change_Style(bpy.types.Operator):
    bl_idname = 'mn.style_change'
    bl_label = 'Style'
    
    style: bpy.props.EnumProperty(
        name = "Style", 
        items = nodes.STYLE_ITEMS
    )
    
    def execute(self, context):
        object = context.active_object
        nodes.change_style_node(object, self.style)
        
        return {'FINISHED'}

def check_installs(selection):
    for package in packages[selection]:
        if not pkg.is_current(package):
            return False
    
    return True

def panel_import(layout, context):
    scene = context.scene
    selection = scene.MN_panel_import
    layout.prop(scene, 'MN_panel_import')
    
    apple_warning = pkg._is_apple_silicon and selection == "md" and not pkg.is_current('MDAnalysis')
    install_required = not check_installs(selection) or apple_warning
    
    if apple_warning:
        pref.apple_silicon_warning(layout)
    
    buttons = layout.column(align=True)
    buttons.enabled = not apple_warning
    
    if install_required:
        buttons.label(text = 'Please install the requried packages.')
        for package in packages[selection]:
            pkg.button_install_pkg(buttons, package, pkg.get_pkgs()[package]['version'])
    
    col = layout.column()
    col.enabled = not install_required
    chosen_panel[selection].panel(col, scene)


def ui_from_node(layout, node):
    """
    Generate the UI for a particular node, which displays the relevant node inputs
    for user control in a panel, rather than through the node editor.
    """
    col = layout.column(align = True)
    ntree = bpy.context.active_object.modifiers['MolecularNodes'].node_group
    
    tree = node.node_tree.interface.items_tree
    
    for item in tree.values():
        if item.item_type == "PANEL":
            col.label(text=item.name)
        elif item.name == "Selection":
            continue
        else:
            if item.in_out != "INPUT":
                continue
            if item.socket_type == "NodeSocketGeometry":
                continue
            # col.prop(node.inputs[item.identifier], 'default_value', text = item.name)
            col.template_node_view(ntree, node, node.inputs[item.identifier])

def panel_object(layout, context):
    scene = context.scene
    object = context.active_object
    mol_type = object.mn.molecule_type
    if mol_type == "":
        layout.label(text = "No MN object selected.")
        return None
    
    if mol_type == "pdb":
        layout.label(text = f"PDB: {object.mn.pdb_code.upper()}")
    if mol_type == "md":
        layout.prop(object.mn, 'subframes')
    if mol_type == "star":
        layout.label(text = f"{object.mn.star_type} starfile")
        box = layout.box()
        ui_from_node(box, nodes.get_star_node(object))
        return
        
    
    row = layout.row(align=True)
    row.label(text = "Style")
    current_style = nodes.format_node_name(nodes.get_style_node(object).node_tree.name).replace("Style ", "")
    row.operator_menu_enum('mn.style_change', 'style', text = current_style)
    box = layout.box()
    ui_from_node(box, nodes.get_style_node(object))
    # layout.label(text='Color')
    # box = layout.box()
    # ui_from_node(box, nodes.get_color_node(object))
    row = layout.row()
    row.label(text="Experimental", icon_value=2)
    row.operator('mn.add_armature')


def panel_scene(layout, context):
    scene = context.scene
    
    cam = bpy.data.cameras[bpy.data.scenes["Scene"].camera.name]
    world_shader = bpy.data.worlds["World Shader"].node_tree.nodes["MN_world_shader"]
    grid = layout.grid_flow()
    col = grid.column()
    col.label(text = "World Settings")
    world = col.box()
    world.prop(bpy.data.scenes["Scene"].render, "engine")
    if scene.render.engine == "CYCLES":
        world.prop(bpy.data.scenes["Scene"].cycles, "samples")
    else:
        world.prop(bpy.data.scenes["Scene"].eevee, "taa_render_samples")
    world.label(text = "Background")
    world.prop(world_shader.inputs[1], 'default_value', text = "World Lighting")
    row = world.row()
    row.prop(scene.render, 'film_transparent')
    row.prop(world_shader.inputs[2], 'default_value', text = "")
    
    col = grid.column()
    col.label(text="Camera Settings")
    camera = col.box()
    camera.prop(cam, "lens")
    col = camera.column(align=True)
    row = col.row(align=True)
    row.prop(bpy.data.scenes["Scene"].render, "resolution_x", text = "X")
    row.prop(bpy.data.scenes["Scene"].render, "resolution_y", text = "Y")
    row = camera.grid_flow()
    row.prop(cam.dof, 'use_dof')
    row.prop(bpy.data.scenes["Scene"].render, "use_motion_blur")
    focus = camera.column()
    focus.enabled = cam.dof.use_dof
    focus.prop(cam.dof, 'focus_object')
    distance = focus.row()
    distance.enabled = (cam.dof.focus_object is None)
    distance.prop(cam.dof, 'focus_distance')
    focus.prop(cam.dof, 'aperture_fstop')


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
        row = layout.row(align=True)
        for p in ['import', 'object', 'scene']:
            row.prop_enum(scene, 'MN_panel', p)
        
        # the possible panel functions to choose between
        which_panel = {
            "import": panel_import,
            "scene": panel_scene,
            "object": panel_object
        }
        # call the required panel function with the layout and context
        which_panel[scene.MN_panel](layout, context)
