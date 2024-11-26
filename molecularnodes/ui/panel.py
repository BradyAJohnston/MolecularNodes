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
STYLE_ITEMS = (
    ("spheres", "Spheres", "Space-filling atoms style."),
    ("cartoon", "Cartoon", "Secondary structure cartoons"),
    ("surface", "Surface", "Solvent-accsible surface."),
    ("ribbon", "Ribbon", "Continuous backbone ribbon."),
    ("sticks", "Sticks", "Sticks for each bond."),
    ("ball_and_stick", "Ball and Stick", "Spheres for atoms, sticks for bonds"),
    ("preset_1", "Preset 1", "A pre-made combination of different styles"),
    ("preset_2", "Preset 2", "A pre-made combination of different styles"),
    ("preset_3", "Preset 3", "A pre-made combination of different styles"),
    ("preset_4", "Preset 4", "A pre-made combination of different styles"),
)

bpy.types.Scene.MN_import_style = bpy.props.EnumProperty(
    name="Style",
    description="Default style for importing molecules.",
    items=STYLE_ITEMS,
    default="spheres",
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
    universe = session.trajectories.get(obj.mn.uuid)
    trajectory_is_linked = bool(universe)
    col = layout.column()
    col.enabled = False
    if not trajectory_is_linked:
        col.enabled = True
        col.label(text="Object not linked to a trajectory, please reload one")
        col.prop(obj.mn, "filepath_topology")
        col.prop(obj.mn, "filepath_trajectory")
        col.operator("mn.reload_trajectory")
        return None

    layout.label(text="Trajectory Playback", icon="OPTIONS")
    row = layout.row()
    col = row.column()
    col.prop(obj.mn, "average")
    col.prop(obj.mn, "subframes")
    col.prop(obj.mn, "offset")
    col = row.column()

    # only enable this as an option if the universe is orthothombic
    row = col.row()
    row.prop(obj.mn, "correct_periodic")
    row.enabled = universe.is_orthorhombic
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
    obj = context.active_object
    session = get_session()
    ent = session.get(obj.mn.uuid)

    layout.label(text=obj.mn.entity_type.title())

    if obj.mn.entity_type == "":
        layout.prop(obj.mn, "entity_type")

    elif obj.mn.entity_type == "trajectory":
        panel_md_properties(layout, context)
    elif obj.mn.entity_type == "molecule":
        layout.label(text=f"Molecule: {obj.name}")
        layout.label(text=f"PDB: {obj.mn.pdb_code}")
    elif obj.mn.entity_type == "ensemble":
        box = layout.box()
        ui_from_node(box, nodes.get_star_node(obj), context=context)
    else:
        layout.label(text="No MN molecular entity selected")


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


class MN_UL_SessionListUI(bpy.types.UIList):
    def filter_items(self, context, data, propname):
        # Get collection property being filtered
        items = getattr(data, propname)

        # Initialize flag array with correct length
        flags = [self.bitflag_filter_item] * len(items)

        # Filter based on entity_type
        for i, item in enumerate(items):
            if item.mn.entity_type == "":
                flags[i] &= ~self.bitflag_filter_item

        # No ordering, return None for order list
        return flags, []

    def draw_item(
        self, context, layout, data, item, icon, active_data, active_propname, index
    ):
        custom_icon = "VIS_SEL_11"

        if self.layout_type in {"DEFAULT", "COMPACT"}:
            if item.mn.molecule_type == "":
                return
            layout.label(text=item.name)

        elif self.layout_type in {"GRID"}:
            layout.alignment = "CENTER"
            layout.label(text="", icon=custom_icon)


def panel_session_list(layout: bpy.types.UILayout, context):
    layout.label(text="testing")
    layout.template_list(
        "MN_UL_SessionListUI",
        "A list",
        bpy.data,
        "objects",
        bpy.context.scene,
        "MN_session_panel_selection",
        rows=5,
    )


def panel_session(layout, context):
    session = get_session(context)
    panel_session_list(layout, context)
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


class MN_Panel_Object(bpy.types.Panel):
    bl_label = "Molecular Nodes"

    def draw(self, context):
        layout = self.layout
        panel_object(layout, context)


class MN_PT_Object(MN_Panel_Object):
    bl_idname = "MN_PT_panel_object"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "object"


class MN_PT_VIEW3D(MN_Panel_Object):
    bl_idname = "MN_PT_panel_view3d"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Item"


CLASSES = [MN_PT_Scene, MN_PT_Object, MN_PT_VIEW3D, MN_UL_SessionListUI]

# ('WINDOW', 'HEADER', 'CHANNELS', 'TEMPORARY', 'UI', 'TOOLS', 'TOOL_PROPS', 'ASSET_SHELF',
# 'ASSET_SHELF_HEADER', 'PREVIEW', 'HUD', 'NAVIGATION_BAR', 'EXECUTE', 'FOOTER',
# 'TOOL_HEADER', 'XR')
