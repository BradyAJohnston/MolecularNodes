import bpy
from ..blender import nodes
from ..io import alphafold, cellpack, density, dna, local, md, star, wwpdb

bpy.types.Scene.MN_panel = bpy.props.EnumProperty(
    name="Panel Selection",
    items=(
        ("import", "Import", "Import macromolecules", 0),
        ("object", "Object", "Adjust settings affecting the selected object", 1),
        ("scene", "Scene", "Change settings for the world and rendering", 2),
    ),
)

bpy.types.Scene.MN_panel_import = bpy.props.EnumProperty(
    name="Method",
    items=(
        ("pdb", "PDB", "Download from the PDB"),
        ("alphafold", "AlphaFold", "Download from the AlphaFold DB"),
        ("local", "Local", "Open a local file"),
        ("md", "MD", "Import a molecular dynamics trajectory"),
        ("density", "Density", "Import an EM Density Map"),
        ("star", "Starfile", "Import a .starfile mapback file"),
        ("cellpack", "CellPack", "Import a CellPack .cif/.bcif file"),
        ("dna", "oxDNA", "Import an oxDNA file"),
    ),
)
STYLE_ITEMS = (
    ("presets", "Presets", "A pre-made combination of different styles"),
    ("spheres", "Spheres", "Space-filling atoms style."),
    ("surface", "Surface", "Solvent-accsible surface."),
    ("cartoon", "Cartoon", "Secondary structure cartoons"),
    ("ribbon", "Ribbon", "Continuous backbone ribbon."),
    ("sticks", "Sticks", "Sticks for each bond."),
    ("ball_and_stick", "Ball and Stick", "Spheres for atoms, sticks for bonds"),
)

bpy.types.Scene.MN_import_style = bpy.props.EnumProperty(
    name="Style",
    description="Default style for importing molecules.",
    items=STYLE_ITEMS,
    default="spheres",
)

chosen_panel = {
    "pdb": wwpdb,
    "local": local,
    "alphafold": alphafold,
    "star": star,
    "md": md,
    "density": density,
    "cellpack": cellpack,
    "dna": dna,
}

packages = {
    "pdb": ["biotite"],
    "alphafold": ["biotite"],
    "star": ["starfile", "mrcfile", "pillow"],
    "local": ["biotite"],
    "cellpack": ["biotite", "msgpack"],
    "md": ["MDAnalysis"],
    "density": ["mrcfile"],
    "dna": [],
}


def change_style_menu(self, context):
    layout = self.layout
    # bob = context.active_object
    layout.label(text="Molecular Nodes")

    # current_style = nodes.get_style_node(bob).replace("Style ", "")
    layout.operator_menu_enum("mn.style_change", "style", text="Style")
    layout.separator()


def is_style_node(context):
    node = context.space_data.edit_tree.nodes.active
    return node.name.startswith("Style")


def change_style_node_menu(self, context):
    layout = self.layout
    layout.label(text="Molecular Nodes", icon="MOD_PARTICLES")
    row = layout.row()
    if is_style_node(context):
        node = context.active_node
        row.operator_menu_enum("mn.style_change_node", "style", text="Change Style")

    # layout.row().column().prop(
    #     context.space_data.edit_tree.nodes.active.node_tree, "color_tag"
    # )

    layout.separator()


def panel_import(layout, context):
    scene = context.scene
    selection = scene.MN_panel_import
    layout.prop(scene, "MN_panel_import")
    buttons = layout.column(align=True)

    col = layout.column()
    chosen_panel[selection].panel(col, scene)


def ui_from_node(layout, node):
    """
    Generate the UI for a particular node, which displays the relevant node inputs
    for user control in a panel, rather than through the node editor.
    """
    col = layout.column(align=True)
    ntree = bpy.context.active_object.modifiers["MolecularNodes"].node_group

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
            col.template_node_view(ntree, node, node.inputs[item.identifier])


def panel_custom_selections(layout, context):
    layout.label(text="Custom Selections")
    scene = context.scene
    bob = context.active_object
    row = layout.row(align=True)

    row = row.split(factor=0.9)
    row.template_list(
        "MN_UL_TrajectorySelectionListUI",
        "A list",
        context.active_object,
        "mn_universe_selections",
        scene,
        "list_index",
        rows=3,
    )
    col = row.column()
    col.operator("mda.new_item", icon="ADD", text="")
    col.operator("mda.delete_item", icon="REMOVE", text="")
    if bob.mn_universe_selections:
        item = bob.mn_universe_selections[bob.mn.universe_selection_index]

        col = layout.column(align=False)
        col.separator()

        row = col.row()
        row.prop(item, "name")
        row.prop(item, "updating")
        row.prop(item, "periodic")
        col.prop(item, "selection")


def panel_object(layout, context):
    object = context.active_object
    mol_type = object.mn.molecule_type
    if mol_type == "":
        layout.label(text="No MN object selected")
        return None

    if mol_type == "pdb":
        layout.label(text=f"PDB: {object.mn.pdb_code.upper()}")
    if mol_type == "md":
        layout.prop(object.mn, "subframes")
        panel_custom_selections(layout, context)
    if mol_type == "star":
        layout.label(text=f"Ensemble")
        box = layout.box()
        ui_from_node(box, nodes.get_star_node(object))
        return

    layout.separator()
    row = layout.row(align=True)
    row.label(text="Style")
    current_style = nodes.get_style_node(object).node_tree.name.replace("Style ", "")
    row.operator_menu_enum("mn.style_change", "style", text=current_style)
    box = layout.box()
    ui_from_node(box, nodes.get_style_node(object))
    row = layout.row()
    # row.label(text="Experimental", icon_value=2)
    # row.operator("mn.add_armature")


def panel_scene(layout, context):
    scene = context.scene

    cam = bpy.data.cameras[bpy.data.scenes["Scene"].camera.name]
    world_shader = bpy.data.worlds["World Shader"].node_tree.nodes["MN_world_shader"]
    grid = layout.grid_flow()
    col = grid.column()
    col.label(text="World Settings")
    world = col.box()
    world.prop(bpy.data.scenes["Scene"].render, "engine")
    if scene.render.engine == "CYCLES":
        world.prop(bpy.data.scenes["Scene"].cycles, "samples")
    else:
        world.prop(bpy.data.scenes["Scene"].eevee, "taa_render_samples")
    world.label(text="Background")
    world.prop(world_shader.inputs[1], "default_value", text="HDRI Strength")
    row = world.row()
    row.prop(scene.render, "film_transparent")
    row.prop(world_shader.inputs[2], "default_value", text="Background")

    col = grid.column()
    col.label(text="Camera Settings")
    camera = col.box()
    camera.prop(cam, "lens")
    col = camera.column(align=True)
    row = col.row(align=True)
    row.prop(bpy.data.scenes["Scene"].render, "resolution_x", text="X")
    row.prop(bpy.data.scenes["Scene"].render, "resolution_y", text="Y")
    row = camera.grid_flow()
    row.prop(cam.dof, "use_dof")
    row.prop(bpy.data.scenes["Scene"].render, "use_motion_blur")
    focus = camera.column()
    focus.enabled = cam.dof.use_dof
    focus.prop(cam.dof, "focus_object")
    distance = focus.row()
    distance.enabled = cam.dof.focus_object is None
    distance.prop(cam.dof, "focus_distance")
    focus.prop(cam.dof, "aperture_fstop")


class MN_PT_panel(bpy.types.Panel):
    bl_label = "Molecular Nodes"
    bl_idname = "MN_PT_panel"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "scene"
    bl_order = 0
    bl_options = {"HEADER_LAYOUT_EXPAND"}
    bl_ui_units_x = 0

    def draw(self, context):
        layout = self.layout
        scene = context.scene
        row = layout.row(align=True)
        for p in ["import", "object", "scene"]:
            row.prop_enum(scene, "MN_panel", p)

        # the possible panel functions to choose between
        which_panel = {
            "import": panel_import,
            "scene": panel_scene,
            "object": panel_object,
        }
        # call the required panel function with the layout and context
        which_panel[scene.MN_panel](layout, context)


CLASSES = []
