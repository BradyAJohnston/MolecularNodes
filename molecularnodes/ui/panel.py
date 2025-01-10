import bpy

from ..entities.trajectory import dna

from ..blender import nodes
from ..session import get_session
from ..entities import density, ensemble, molecule, trajectory

bpy.types.Scene.MN_panel = bpy.props.EnumProperty(
    name="Panel Selection",
    items=(
        ("import", "Import", "Import macromolecules", 0),
        ("object", "Object", "Adjust settings affecting the selected object", 1),
        (
            "session",
            "Session",
            "Interacting with the Molecular Nodes session tracking all of the objects",
            2,
        ),
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


chosen_panel = {
    "pdb": molecule.ui.panel_wwpdb,
    "local": molecule.ui.panel_local,
    "alphafold": molecule.ui.panel_alphafold,
    "star": ensemble.ui.panel_starfile,
    "md": trajectory.ui.panel,
    "density": density.ui.panel,
    "cellpack": ensemble.ui.panel_cellpack,
    "dna": dna.panel,
}


def pt_object_context(self, context):
    layout = self.layout
    return None


def is_style_node(context):
    node = context.space_data.edit_tree.nodes.active
    return node.name.startswith("Style")


def change_style_node_menu(self, context):
    layout = self.layout
    node = context.active_node
    prefix = node.node_tree.name.split(" ")[0].lower()
    if prefix not in ["color", "select", "is", "style", "topology", "animate"]:
        return None
    layout.label(text="Molecular Nodes", icon="MOD_PARTICLES")

    row = layout.row()
    op = row.operator_menu_enum("mn.node_swap", "node_items", text="Change Node")
    op.node_description = "The topology nodes"

    layout.separator()


def panel_import(layout, context):
    scene = context.scene
    selection = scene.MN_panel_import
    layout.prop(scene, "MN_panel_import")

    col = layout.column()
    chosen_panel[selection](col, scene)


def ui_from_node(
    layout: bpy.types.UILayout, node: bpy.types.NodeGroup, context: bpy.types.Context
):
    """
    Generate the UI for a particular node, which displays the relevant node inputs
    for user control in a panel, rather than through the node editor.
    """
    col = layout.column(align=True)
    ntree = context.active_object.modifiers["MolecularNodes"].node_group

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


def panel_md_properties(layout, context):
    obj = context.active_object
    session = get_session()
    traj: trajectory.Trajectory = session.match(obj)
    traj_is_linked = bool(traj)
    if traj is not None and not isinstance(traj, trajectory.Trajectory):
        raise TypeError(f"Expected a trajectory, got {type(traj)}")

    col = layout.column()
    col.enabled = False
    if not traj_is_linked:
        col.enabled = True
        col.label(text="Object not linked to a trajectory, please reload one")
        col.prop(obj.mn, "filepath_topology")
        col.prop(obj.mn, "filepath_trajectory")
        col.operator("mn.reload_trajectory")
        return None

    layout.label(text="Trajectory Playback", icon="OPTIONS")

    row = layout.row()
    col = row.column()
    if obj.mn.update_with_scene:
        col.prop(obj.mn, "frame_hidden")
    else:
        col.prop(obj.mn, "frame")
    col.enabled = not obj.mn.update_with_scene
    row.prop(obj.mn, "update_with_scene")
    row = layout.row()
    col = row.column()
    col.enabled = obj.mn.update_with_scene
    col.prop(obj.mn, "average")
    col.prop(obj.mn, "subframes")
    col.prop(obj.mn, "offset")
    col = row.column()
    col.enabled = obj.mn.update_with_scene

    # only enable this as an option if the universe is orthothombic
    row = col.row()
    row.prop(obj.mn, "correct_periodic")
    row.enabled = traj.is_orthorhombic
    col.prop(obj.mn, "interpolate")

    layout.label(text="Selections", icon="RESTRICT_SELECT_OFF")
    row = layout.row()
    row = row.split(factor=0.9)
    row.template_list(
        "MN_UL_TrajectorySelectionListUI",
        "A list",
        obj,
        "mn_trajectory_selections",
        obj.mn,
        "trajectory_selection_index",
        rows=3,
    )
    col = row.column()
    col.operator("mn.trajectory_selection_add", icon="ADD", text="")
    col.operator("mda.delete_item", icon="REMOVE", text="")
    if obj.mn_trajectory_selections:
        item = obj.mn_trajectory_selections[obj.mn.trajectory_selection_index]

        col = layout.column(align=False)
        row = col.row()
        col.prop(item, "selection_str")

        # disable editing for immutable selections
        # disable modifying updating and periodic
        if item.immutable:
            col.enabled = False

        if item.message != "":
            box = col.box()
            box.label(text="Invalid Selection", icon="ERROR")
            box.label(text=item.message)
            box.alert = True
            op = box.operator("wm.url_open", text="Selection Langauge Docs", icon="URL")
            op.url = (
                "https://docs.mdanalysis.org/stable/documentation_pages/selections.html"
            )


def panel_object(layout, context):
    object = context.active_object
    layout.prop(object.mn, "entity_type")
    try:
        mol_type = object.mn.entity_type
    except AttributeError:
        return None
    if mol_type == "":
        layout.label(text="No MN object selected")
        return None
    if mol_type == "pdb":
        layout.label(text=f"PDB: {object.mn.pdb_code.upper()}")
    if mol_type.startswith("md"):
        panel_md_properties(layout, context)
    if mol_type == "star":
        layout.label(text="Ensemble")
        box = layout.box()
        ui_from_node(box, nodes.get_star_node(object), context=context)
        return None


def item_ui(layout, item):
    row = layout.row()
    row.label(text=item.name)
    col = row.column()
    op = col.operator("mn.session_create_object")
    op.uuid = item.uuid
    col.enabled = item.object is None

    op = row.operator("mn.session_remove_item", text="", icon="CANCEL")
    op.uuid = item.uuid

    if item.object is not None:
        row = layout.row()
        row.label(text=f"Object: {item.object.name}", icon="OUTLINER_OB_MESH")


def panel_session(layout, context):
    session = get_session(context)
    # if session.n_items > 0:
    #     return None
    row = layout.row()
    row.label(text="Loaded items in the session")
    # row.operator("mn.session_reload")

    layout.label(text="Molecules")
    box = layout.box()
    for mol in session.molecules.values():
        item_ui(box, mol)

    layout.label(text="Universes")
    box = layout.box()
    for uni in session.trajectories.values():
        item_ui(box, uni)

    layout.label(text="Ensembles")
    box = layout.box()
    for ens in session.ensembles.values():
        item_ui(box, ens)


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


class MN_PT_Scene(bpy.types.Panel):
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
        for p in ["import", "object", "session"]:
            row.prop_enum(scene, "MN_panel", p)

        # the possible panel functions to choose between
        which_panel = {
            "import": panel_import,
            "object": panel_object,
            "session": panel_session,
        }
        # call the required panel function with the layout and context
        which_panel[scene.MN_panel](layout, context)


CLASSES = []
