import bpy
from ..entities import trajectory
from ..nodes import nodes
from ..session import get_session
from .pref import addon_preferences
from .utils import check_online_access_for_ui


def panel_wwpdb(layout, scene):
    layout.label(text="Download from PDB", icon="IMPORT")
    layout.separator()

    layout = check_online_access_for_ui(layout)

    row_import = layout.row().split(factor=0.5)
    row_import.prop(scene.mn, "import_code_pdb")
    row = row_import.split(factor=0.3)
    row.prop(scene.mn, "import_format_wwpdb", text="")
    op = row.operator("mn.import_fetch")
    op.code = scene.mn.import_code_pdb  # type: ignore
    op.database = "wwpdb"  # type: ignore
    op.file_format = scene.mn.import_format_wwpdb  # type: ignore
    op.node_setup = scene.mn.import_node_setup  # type: ignore
    op.remove_solvent = scene.mn.import_remove_solvent  # type: ignore
    op.assembly = scene.mn.import_build_assembly  # type: ignore
    op.style = scene.mn.import_style  # type: ignore
    op.centre = scene.mn.import_centre  # type: ignore
    op.centre_type = scene.mn.import_centre_type  # type: ignore
    prefs = addon_preferences()
    if prefs is not None:
        op.cache_dir = str(prefs.cache_dir)  # type: ignore
    else:
        op.cache_dir = str(bpy.app.tempdir)  # type: ignore
    layout.separator(factor=0.4)

    layout.separator()

    layout.label(text="Options", icon="MODIFIER")
    options = layout.column(align=True)

    row = options.row()
    row.prop(scene.mn, "import_node_setup", text="")
    col = row.column()
    col.prop(scene.mn, "import_style")
    col.enabled = scene.mn.import_node_setup

    row_centre = options.row()
    row_centre.prop(scene.mn, "import_centre", icon_value=0)
    col_centre = row_centre.column()
    col_centre.prop(scene.mn, "import_centre_type", text="")
    col_centre.enabled = scene.mn.import_centre
    options.separator()

    grid = options.grid_flow()
    grid.prop(scene.mn, "import_build_assembly")
    grid.prop(scene.mn, "import_remove_solvent")
    grid.prop(scene.mn, "import_del_hydrogen")


def panel_alphafold(layout, scene):
    layout.label(text="Download from the AlphaFold DataBase", icon="IMPORT")
    layout.separator()

    layout = check_online_access_for_ui(layout)

    row_import = layout.row().split(factor=0.5)
    row_import.prop(scene.mn, "import_code_alphafold")
    download = row_import.split(factor=0.3)
    download.prop(scene.mn, "import_format_alphafold", text="")
    op = download.operator("mn.import_fetch")
    op.code = scene.mn.import_code_alphafold  # type: ignore
    op.database = "alphafold"  # type: ignore
    op.file_format = scene.mn.import_format_alphafold  # type: ignore
    op.node_setup = scene.mn.import_node_setup  # type: ignore
    op.assembly = scene.mn.import_build_assembly  # type: ignore
    op.style = scene.mn.import_style  # type: ignore
    op.centre = scene.mn.import_centre  # type: ignore
    op.centre_type = scene.mn.import_centre_type  # type: ignore
    prefs = addon_preferences()
    if prefs is not None:
        op.cache_dir = str(prefs.cache_dir)  # type: ignore
    else:
        op.cache_dir = str(bpy.app.tempdir)  # type: ignore

    layout.separator(factor=0.4)

    row = layout.row().split(factor=0.3)
    layout.separator()

    layout.label(text="Options", icon="MODIFIER")
    options = layout.column(align=True)

    row = options.row()
    row.prop(scene.mn, "import_node_setup", text="")
    col = row.column()
    col.prop(scene.mn, "import_style")
    col.enabled = scene.mn.import_node_setup

    row_centre = options.row()
    row_centre.prop(scene.mn, "import_centre", icon_value=0)
    col_centre = row_centre.column()
    col_centre.prop(scene.mn, "centre_type", text="")
    col_centre.enabled = scene.mn.import_centre
    options.separator()


# operator that calls the function to import the structure from a local file


def panel_local(layout, scene):
    layout.label(text="Load a Local File", icon="FILE_TICK")
    layout.separator()

    row = layout.row()
    row.prop(scene.mn, "import_local_path")
    op = row.operator("mn.import_local")
    op.filepath = scene.mn.import_local_path
    op.node_setup = scene.mn.import_node_setup
    op.assembly = scene.mn.import_build_assembly
    op.style = scene.mn.import_style
    op.centre = scene.mn.import_centre
    op.remove_solvent = scene.mn.import_remove_solvent
    op.centre_type = scene.mn.import_centre_type
    layout.separator()

    layout.label(text="Options", icon="MODIFIER")
    options = layout.column(align=True)

    row = options.row()
    row.prop(scene.mn, "import_node_setup", text="")
    col = row.column()
    col.prop(scene.mn, "import_style")
    col.enabled = scene.mn.import_node_setup

    row_centre = options.row()

    row_centre.prop(scene.mn, "import_centre", icon_value=0)
    # row_centre.prop()
    col_centre = row_centre.column()
    col_centre.prop(scene.mn, "import_centre_type", text="")
    col_centre.enabled = scene.mn.import_centre
    options.separator()

    grid = options.grid_flow()
    grid.prop(scene.mn, "import_build_assembly")
    grid.prop(scene.mn, "import_remove_solvent", icon_value=0)
    grid.prop(scene.mn, "import_del_hydrogen", icon_value=0)


def panel_starfile(layout, scene):
    layout.label(text="Load Star File", icon="FILE_TICK")
    layout.separator()
    row_import = layout.row()
    row_import.prop(scene.mn, "import_star_file_path")
    op = row_import.operator("mn.import_star_file")
    op.filepath = scene.mn.import_star_file_path
    op.node_setup = scene.mn.import_node_setup


def panel_cellpack(layout, scene):
    layout.label(text="Load CellPack Model", icon="FILE_TICK")
    layout.separator()
    row = layout.row()
    row.prop(scene.mn, "import_cell_pack_path")
    op = row.operator("mn.import_cell_pack")
    op.filepath = scene.mn.import_cell_pack_path
    op.node_setup = scene.mn.import_node_setup


def panel_density(layout, scene):
    layout.label(text="Load EM Map", icon="FILE_TICK")
    layout.separator()

    row = layout.row()
    row.prop(scene.mn, "import_density")
    row.operator("mn.import_density")

    layout.separator()
    col = layout.column()
    col.alignment = "LEFT"
    col.scale_y = 0.5
    label = f"\
    An intermediate file will be created: {scene.mn.import_density}.vdb\
    Please do not delete this file or the volume will not render.\
    Move the original .map file to change this location.\
    "
    for line in label.strip().split("    "):
        col.label(text=line)

    layout.separator()
    layout.label(text="Options", icon="MODIFIER")

    layout.prop(scene.mn, "import_density_invert")
    layout.prop(scene.mn, "import_density_center")
    row = layout.row()
    row.prop(scene.mn, "import_node_setup", text="")
    col = row.column()
    col.prop(scene.mn, "import_density_style")
    col.enabled = scene.mn.import_node_setup


def panel_trajectory(layout, scene):
    layout.label(text="Load MD Trajectories", icon="FILE_TICK")
    layout.separator()
    col = layout.column(align=True)
    row_import = col.row()
    row_import.prop(scene.mn, "import_md_name")
    op = row_import.operator("mn.import_trajectory", text="Load")
    op.topology = scene.mn.import_md_topology
    op.trajectory = scene.mn.import_md_trajectory
    op.name = scene.mn.import_md_name
    op.style = scene.mn.import_style
    op.setup_nodes = scene.mn.import_node_setup
    col.separator()
    col.prop(scene.mn, "import_md_topology")
    col.prop(scene.mn, "import_md_trajectory")

    layout.separator()
    layout.label(text="Options", icon="MODIFIER")
    row = layout.row()
    row.prop(scene.mn, "import_node_setup", text="")
    col = row.column()
    col.prop(scene.mn, "import_style")
    col.enabled = scene.mn.import_node_setup


def panel_oxdna(layout: bpy.types.UILayout, scene: bpy.types.Scene) -> None:
    """
    Create the panel layout for oxDNA import.

    Parameters
    ----------
    layout : bpy.types.UILayout
        Layout to add elements to
    scene : bpy.types.Scene
        Current scene
    """
    layout.label(text="Load oxDNA File", icon="FILE_TICK")
    layout.separator()
    row = layout.row()
    row.prop(scene.mn, "import_oxdna_name")
    op = row.operator("mn.import_oxdna")
    op.name = scene.mn.import_oxdna_name
    op.topology = scene.mn.import_oxdna_topology
    op.trajectory = scene.mn.import_oxdna_trajectory
    col = layout.column(align=True)
    col.prop(scene.mn, "import_oxdna_topology")
    col.prop(scene.mn, "import_oxdna_trajectory")


chosen_panel = {
    "pdb": panel_wwpdb,
    "local": panel_local,
    "alphafold": panel_alphafold,
    "star": panel_starfile,
    "md": panel_trajectory,
    "density": panel_density,
    "cellpack": panel_cellpack,
    "dna": panel_oxdna,
}


def pt_object_context(self, context):
    return None


def is_style_node(context):
    node = context.space_data.edit_tree.nodes.active
    return node.name.startswith("Style")


def change_style_node_menu(self, context):
    layout = self.layout
    node = context.active_node

    # return early if not a node group
    if not hasattr(node, "node_tree"):
        return None

    # return early if the node group isn't one of the ones we want to swap easily
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
    selection = scene.mn.panel_import_type
    layout.prop(scene.mn, "panel_import_type")

    col = layout.column()
    chosen_panel[selection](col, scene)


def ui_from_node(
    layout: bpy.types.UILayout,
    node: bpy.types.GeometryNodeGroup,
    context: bpy.types.Context,
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
    label = "This trajectory has " + str(obj.mn.n_frames) + " frame"
    if obj.mn.n_frames > 1:
        label += "s"
    layout.label(text=label)

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
    if object is None:
        # When an object is deleted, context.ative_object is None
        return
    layout.prop(object.mn, "entity_type")
    try:
        mol_type = object.mn.entity_type
    except AttributeError:
        return None
    if mol_type == "":
        layout.label(text="No MN object selected")
        return None
    if mol_type == "pdb":
        layout.label(text=f"PDB: {object.mn.code.upper()}")
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
    row = layout.row()
    row.label(text="Loaded items in the session")

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
        row = layout.row()

        row = layout.row(align=True)

        for p in ["import", "object", "session"]:
            row.prop_enum(scene.mn, "panel_selection", p)

        # the possible panel functions to choose between
        which_panel = {
            "import": panel_import,
            "object": panel_object,
            "session": panel_session,
        }
        # call the required panel function with the layout and context
        which_panel[scene.mn.panel_selection](layout, context)


CLASSES = [MN_PT_Scene]
