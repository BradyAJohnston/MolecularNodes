from typing import cast
import bpy
from bpy.types import UILayout
from ..entities import Molecule, StreamingTrajectory, molecule
from ..entities.base import EntityType
from ..nodes import nodes
from ..session import get_session
from .ops import MN_OT_add_selection_to_style
from .props import TrajectorySelectionItem
from .utils import check_online_access_for_ui


def import_options(layout: UILayout) -> None:
    # fetching requires online access, so gate only those entries
    online = check_online_access_for_ui(layout.column())
    op = online.operator("mn.import_molecule", text="Fetch from PDB", icon="IMPORT")
    op.method = "fetch"
    op.database = "wwpdb"
    op = online.operator(
        "mn.import_molecule", text="Fetch from AlphaFold", icon="IMPORT"
    )
    op.database = "alphafold"

    op = layout.operator("mn.import_molecule", text="Import Local File", icon="IMPORT")
    op.method = "local"
    op = layout.operator("mn.import_molecule", text="From SMILES", icon="IMPORT")
    op.method = "smiles"

    layout.separator()
    op = layout.operator("mn.import_ensemble", text="Starfile", icon="IMPORT")
    op.ensemble_type = "starfile"
    op = layout.operator("mn.import_ensemble", text="CellPack", icon="IMPORT")
    op.ensemble_type = "cellpack"

    layout.separator()
    # MD trajectories are loaded through the same unified import dialog
    op = layout.operator("mn.import_molecule", text="MD Trajectory", icon="IMPORT")
    op.method = "local"
    op = layout.operator("mn.import_oxdna", text="oxDNA", icon="IMPORT")

    layout.separator()
    layout.operator("mn.import_density", text="Density Map", icon="IMPORT")


class MN_MT_Add(bpy.types.Menu):
    bl_idname = "MN_MT_Add"
    bl_label = "Molecular Nodes"

    def draw(self, context: bpy.types.Context) -> None:
        layout = self.layout
        assert layout
        # the Add menu defaults to EXEC_REGION_WIN, which would run the operator
        # directly; force INVOKE so the operators' popup dialogs are shown instead
        with context.temp_override(operator_context="INVOKE_DEFAULT"):
            import_options(layout)


def add_menu_options(self: bpy.types.Menu, context: bpy.types.Context) -> None:
    layout = self.layout
    assert layout
    layout.menu("MN_MT_Add", icon="PARTICLES")


class MN_MT_Import(bpy.types.Menu):
    bl_idname = "MN_MT_Import"
    bl_label = "Import"

    def draw(self, context: bpy.types.Context) -> None:
        layout = self.layout
        assert layout
        # force INVOKE so the operators' popup dialogs are shown instead of
        # being executed directly
        with context.temp_override(operator_context="INVOKE_DEFAULT"):
            import_options(layout)


def pt_object_context(self, context):
    # currently does and returns nothing but I've left it in so we can
    # know where to add stuff in later for context-dependent menus
    return None


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
    ntree = context.active_object.modifiers["Molecular Nodes"].node_group

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


def selection_string_input(layout: bpy.types.UILayout, item: TrajectorySelectionItem):
    row = layout.row(align=True)
    row.prop(item, "string")

    # disable editing for immutable selections
    # disable modifying updating and periodic
    if item.from_atomgroup:
        row.enabled = False

    if item.message != "":
        box = layout.box()
        box.label(text="Invalid Selection", icon="ERROR")
        box.label(text=item.message)
        box.alert = True
        op = box.operator("wm.url_open", text="Selection Langauge Docs", icon="URL")
        op.url = (
            "https://docs.mdanalysis.org/stable/documentation_pages/selections.html"
        )

    return row


def layout_trajectory_playback(
    layout: UILayout, traj: Molecule, panel: bool = True
) -> None:
    if panel:
        header, layout = layout.panel(idname="layout_playback")
        header.label(text="Trajectory Playback", icon="OPTIONS")
        if layout is None:
            return
    obj = traj.object
    is_streaming = isinstance(traj, StreamingTrajectory)

    # a single static structure (one topology, no trajectory) has nothing to play back,
    # so the playback controls are only shown when there are multiple frames
    if not is_streaming and obj.mn.n_frames <= 1:
        if panel:
            layout.label(text="Single frame; no playback")
        return

    if is_streaming:
        label = "Streaming trajectory; cannot alter playback"
    else:
        label = "Trajectory has {} {}".format(
            obj.mn.n_frames, "frames" if obj.mn.n_frames else "frame"
        )

    layout.label(text=label)
    playback = layout.column()
    playback.enabled = not is_streaming
    row = playback.row()

    col = row.column()
    if obj.mn.update_with_scene:
        col.prop(obj.mn, "frame_hidden")
    else:
        col.prop(obj.mn, "frame")
    col.enabled = not obj.mn.update_with_scene
    row.prop(obj.mn, "update_with_scene")
    row = playback.row()
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
    row.enabled = traj._is_orthorhombic
    col.prop(obj.mn, "interpolate")

    # bake frames into a collection for use directly in geometry nodes
    playback.separator()
    playback.operator("mn.frames_to_collection", icon="RENDERLAYERS")


def layout_selection_manage(
    layout: UILayout, traj: Molecule, panel: bool = True
) -> None:
    if panel:
        header, layout = layout.panel(idname="selection_panel")
        header.label(text="Selections", icon="RESTRICT_SELECT_OFF")
        if layout is None:
            return
    obj = traj.object
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
    col.operator("mn.trajectory_selection_remove", icon="REMOVE", text="")
    if obj.mn_trajectory_selections:
        item = obj.mn_trajectory_selections[obj.mn.trajectory_selection_index]

        selection_string_input(layout, item)


def panel_md_properties(layout, context):
    obj = context.active_object
    session = get_session()
    traj: molecule.Molecule = session.match(obj)
    traj_is_linked = bool(traj)
    if traj is not None and not isinstance(traj, molecule.Molecule):
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

    layout_trajectory_playback(layout, traj)
    layout_selection_manage(layout, traj)


def panel_object(layout: bpy.types.UILayout, context: bpy.types.Context):
    object = cast(bpy.types.Object, context.active_object)
    if object is None:
        # When an object is deleted, context.ative_object is None
        return
    row = layout.row()
    row.prop(object.mn, "entity_type")
    row.enabled = False
    try:
        mol_type = object.mn.entity_type
    except AttributeError:
        return None
    if not object.mn.is_entity:
        layout.label(text="No MN object selected")
        return None
    if mol_type.startswith("md") or mol_type == EntityType.MOLECULE.value:
        # molecules and trajectories are both Universe-backed, so both expose selections
        # and (frame-count permitting) playback
        panel_md_properties(layout, context)
    if mol_type == "ensemble-star":
        layout.label(text="Ensemble")
        box = layout.box()
        ui_from_node(box, nodes.get_star_node(object), context=context)
        return None


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
        assert layout

        # import operators live in a single drop-down menu at the top
        layout.menu("MN_MT_Import", text="Import", icon="IMPORT")

        # display the information for the selected object
        panel_object(layout, context)


class MN_PT_Object(bpy.types.Panel):
    bl_label = "Molecular Nodes"
    bl_idname = "MN_PT_object"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "object"
    bl_order = 0
    bl_options = {"HEADER_LAYOUT_EXPAND"}
    bl_ui_units_x = 0

    def draw(self, context):
        layout = self.layout
        assert layout

        # display the information for the selected object
        panel_object(layout, context)


def get_active_entity_object(context: bpy.types.Context) -> bpy.types.Object | None:
    """
    The object currently selected in the Entities list, or None.

    The Entities list displays `bpy.data.objects` filtered to molecular entities
    (objects with `mn.entity_type` set), so the active index is an index into
    `bpy.data.objects` and can point at a non-entity object after filtering.
    """
    index = context.scene.mn.entities_active_index
    objects = bpy.data.objects
    if 0 <= index < len(objects):
        obj = objects[index]
        if obj.mn.is_entity:
            return obj
    return None


def get_entity_node_group(
    obj: bpy.types.Object | None,
) -> bpy.types.GeometryNodeTree | None:
    """The node tree of the object's "Molecular Nodes" modifier, or None."""
    if obj is None:
        return None
    mod = obj.modifiers.get("Molecular Nodes")
    if mod is None:
        return None
    return getattr(mod, "node_group", None)


class MN_UL_EntitiesList(bpy.types.UIList):
    """
    UIList of molecular entity objects in the Entities panel (Viewport).

    Lists `bpy.data.objects`, filtered down to objects that have
    `mn.entity_type` set. Each object may or may not have a matching entity
    tracked in the session's MNSession, indicated by the link icon.
    """

    def draw_item(
        self,
        context,
        layout,
        data,
        item,
        icon,
        active_data,
        active_property,
        index=0,
        flt_flag=0,
    ):
        obj = cast(bpy.types.Object, item)
        if self.layout_type in {"DEFAULT", "COMPACT"}:
            row = layout.row()
            linked = context.scene.MNSession.get(obj.uuid) is not None
            row.label(text="", icon="LINKED" if linked else "UNLINKED")
            row.prop(obj, "name", text="", emboss=False)
            # use the object viewport visibility to determine the icon
            # we do not have direct callbacks for raw object visibility changes
            hide_icon = "HIDE_OFF" if obj.mn.visible else "HIDE_ON"
            row.prop(obj.mn, "visible", icon_only=True, icon=hide_icon)
        elif self.layout_type in {"GRID"}:
            layout.alignment = "CENTER"
            layout.label(text="", icon="OBJECT_DATA")

    def filter_items(self, context, data, propname):
        if data is None:
            return [], []
        objects = getattr(data, propname)
        # show only molecular entities, then apply the name filter on top
        sort_data = []
        filtered = [0] * len(objects)
        for i, obj in enumerate(objects):
            sort_data.append((i, obj.name))
            if not obj.mn.is_entity:
                continue
            if (
                not self.filter_name
                or (self.filter_name.lower() in obj.name.lower())
                is not self.use_filter_invert
            ):
                filtered[i] |= self.bitflag_filter_item
        # Sort
        ordered = []
        if self.use_filter_sort_alpha:
            ordered = bpy.types.UI_UL_list.sort_items_helper(sort_data, lambda e: e[1])
        return filtered, ordered


class MN_PT_Entities(bpy.types.Panel):
    """
    Panel to list MN Entities in Viewport
    """

    bl_idname = "MN_PT_Entities"
    bl_label = "Entities"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Molecular Nodes"

    def draw(self, context):
        layout = self.layout
        assert layout
        props = context.scene.mn
        row = layout.row()
        row.template_list(
            "MN_UL_EntitiesList",
            "entities_list",
            bpy.data,
            "objects",
            props,
            "entities_active_index",
            rows=3,
        )

        obj = get_active_entity_object(context)
        entity = context.scene.MNSession.get(obj.uuid) if obj is not None else None

        col = row.column()
        # reload/relink the entity into the session from its recorded source
        row = col.row()
        row.operator("mn.session_reload_item", icon="FILE_REFRESH", text="")
        row = col.row()
        op = row.operator("mn.session_remove_item", icon="REMOVE", text="")
        if entity is None:
            row.enabled = False
        else:
            op.uuid = obj.uuid

        if obj is None:
            return
        # display details of the selected entity object
        row = layout.row()
        row.prop(obj.mn, "entity_type")
        row.enabled = False
        row = layout.row()
        if entity is not None:
            row.label(text="Linked to session entity", icon="LINKED")
        else:
            row.label(text="Not linked — use Reload to relink", icon="UNLINKED")


class MN_PT_trajectory(bpy.types.Panel):
    """
    Panel for trajectory details
    """

    bl_idname = "MN_PT_trajectory"
    bl_label = "Trajectory"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Molecular Nodes"

    @classmethod
    def poll(cls, context):
        """Visible for a Universe-backed entity that has something to play back."""
        obj = get_active_entity_object(context)
        if obj is None or context.scene.MNSession.get(obj.uuid) is None:
            return False
        if obj.mn.entity_type not in (
            EntityType.MD.value,
            EntityType.MD_STREAMING.value,
            EntityType.MD_OXDNA.value,
            EntityType.MOLECULE.value,
        ):
            return False
        # a single static structure has no playback; streaming has an unknown frame count
        return (
            obj.mn.entity_type == EntityType.MD_STREAMING.value or obj.mn.n_frames > 1
        )

    def draw(self, context):
        layout = cast(UILayout, self.layout)
        assert layout
        # To enable the animatate dot next to property in UI
        # layout.use_property_split = True
        # layout.use_property_decorate = True
        obj = get_active_entity_object(context)
        traj = context.scene.MNSession.get(obj.uuid)

        layout_trajectory_playback(layout, traj, panel=False)


class MN_PT_trajectory_dssp(bpy.types.Panel):
    """
    Panel for trajectory dssp details
    """

    bl_idname = "MN_PT_trajectory_dssp"
    bl_label = "DSSP"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Molecular Nodes"

    @classmethod
    def poll(cls, context):
        """Visible only if entity selected is a trajectory"""
        obj = get_active_entity_object(context)
        if obj is None or context.scene.MNSession.get(obj.uuid) is None:
            return False
        return obj.mn.entity_type in (
            EntityType.MD.value,
            EntityType.MD_STREAMING.value,
        )

    def draw(self, context):
        layout = self.layout
        assert layout
        obj = get_active_entity_object(context)
        uuid = obj.uuid
        traj = context.scene.MNSession.get(uuid)
        if traj.dssp._DSSP is None:
            row = layout.row()
            op = row.operator("mn.dssp_init")
            op.uuid = uuid
            return
        props = traj.props.dssp
        # display options
        if traj._entity_type == EntityType.MD:
            row = layout.row()
            row.prop(props, "display_option")
        elif traj._entity_type == EntityType.MD_STREAMING:
            row = layout.row()
            row.prop(props, "display_option_streaming")
        # display option specific params
        if props.display_option == "sliding-window-average":
            row = layout.row()
            row.prop(props, "window_size")
            row = layout.row()
            row.prop(props, "apply_sw_threshold", text="")
            col = row.column()
            col.prop(props, "sw_threshold")
            col.enabled = props.apply_sw_threshold
        elif props.display_option == "trajectory-average":
            row = layout.row()
            row.prop(props, "apply_ta_threshold", text="")
            col = row.column()
            col.prop(props, "ta_threshold")
            col.enabled = props.apply_ta_threshold
        # apply button
        if props.display_option == "trajectory-average":
            row = layout.row()
            split = row.split(factor=0.5)
            col = split.column()
            op = col.operator("mn.dssp_apply")
            op.uuid = uuid
            op.apply_ta_threshold = props.apply_ta_threshold
            op.ta_threshold = props.ta_threshold
            col.enabled = not props.applied
            col = split.column()
            op = col.operator("mn.dssp_cancel")
            op.uuid = uuid


def is_style_node(node: bpy.types.Node) -> bool:
    """Whether a node is a style node, identified by "Style" in its name."""
    return "Style" in node.name


class MN_UL_StylesList(bpy.types.UIList):
    """
    UIList of styles for an entity.

    Displays the style nodes of the entity's node tree, filtered directly to the
    nodes with "Style" in their name.
    """

    def draw_item(
        self,
        context,
        layout: bpy.types.UILayout,
        data,
        item,
        icon,
        active_data,
        active_property,
        index=0,
        flt_flag=0,
    ):
        item: bpy.types.GeometryNode = item
        layout: bpy.types.UILayout = layout
        if self.layout_type in {"DEFAULT", "COMPACT"}:
            row = layout.row()
            row.prop(item, "name", text="", emboss=False)
            row.prop(
                item,
                "mute",
                emboss=True,
                icon_only=True,
                icon="RESTRICT_RENDER_OFF" if not item.mute else "RESTRICT_RENDER_ON",
            )
            if "Visible" in item.inputs:
                input = item.inputs["Visible"]
                hide_icon = "HIDE_OFF" if input.default_value else "HIDE_ON"
                row.prop(input, "default_value", icon_only=True, icon=hide_icon)
        elif self.layout_type in {"GRID"}:
            layout.alignment = "CENTER"
            layout.label(text="", icon="WORLD")

    def filter_items(self, context, data, propname):
        if data is None:
            return [], []
        items = getattr(data, propname)
        # show only style nodes, then apply the name filter on top
        sort_data = []
        filtered = [0] * len(items)
        for i, item in enumerate(items):
            if not is_style_node(item):
                sort_data.append((i, ""))
                continue
            name = item.label
            sort_data.append((i, name))
            if (
                not self.filter_name
                or (self.filter_name.lower() in name.lower())
                is not self.use_filter_invert
            ):
                filtered[i] |= self.bitflag_filter_item
        # Sort
        ordered = []
        if self.use_filter_sort_alpha:
            ordered = bpy.types.UI_UL_list.sort_items_helper(sort_data, lambda e: e[1])
        return filtered, ordered


class MN_PT_Styles(bpy.types.Panel):
    """
    Panel for styles
    """

    bl_idname = "MN_PT_styles"
    bl_label = "Styles"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Molecular Nodes"

    @classmethod
    def poll(cls, context):
        """Visible only if the active entity is a trajectory, molecule or density"""
        obj = get_active_entity_object(context)
        return obj is not None and obj.mn.entity_type in (
            EntityType.MD.value,
            EntityType.MD_STREAMING.value,
            EntityType.MOLECULE.value,
            EntityType.DENSITY.value,
        )

    def draw(self, context):
        layout = self.layout
        assert layout is not None
        obj = get_active_entity_object(context)
        if obj is None:
            return
        # style nodes live in the object's node tree — read them directly from
        # Blender data so the panel does not depend on the session being linked
        node_group = get_entity_node_group(obj)
        if node_group is None:
            layout.label(text="No Molecular Nodes modifier on this object")
            return

        # list the style nodes in the tree and let the user select one
        layout.template_list(
            "MN_UL_StylesList",
            "styles_list",
            node_group,
            "nodes",
            obj.mn,
            "styles_active_index",
            rows=3,
        )

        # the style node selected in the list, if any
        index = obj.mn.styles_active_index
        style_node = None
        if 0 <= index < len(node_group.nodes):
            node = node_group.nodes[index]
            if is_style_node(node):
                style_node = node
        if style_node is None:
            layout.label(text="Select a style to edit its properties")
            return

        # swap the selected style node for a different style
        row = layout.row(align=True)
        row.label(text="Style:")
        op = row.operator_menu_enum(
            "mn.swap_style",
            "style",
            text=style_node.node_tree.name.replace("Style ", ""),
        )
        op.name_tree = node_group.name
        op.name_node = style_node.name

        # display the selection string if using a named attribute
        entity = context.scene.MNSession.get(obj.uuid)

        if style_node.inputs["Selection"].links and entity is not None:
            node = style_node.inputs["Selection"].links[0].from_node
            if isinstance(node, bpy.types.GeometryNodeInputNamedAttribute):
                if isinstance(entity, Molecule):
                    selection = entity.selections.get(node.inputs["Name"].default_value)
                    layout.prop(selection, "string", text="Selection")
        else:
            op: MN_OT_add_selection_to_style = layout.operator(
                operator="mn.add_selection_to_style"
            )
            op.node_tree = node_group.name
            op.node_name = style_node.name

        # display the selected style node's name and its input properties
        header, body = layout.panel(idname="style_properties")
        header.label(text=style_node.label or style_node.name)
        if body is not None:
            body.template_node_inputs(style_node)


class MN_UL_AnnotationsList(bpy.types.UIList):
    """
    UIList of annotations for an entity
    """

    def draw_item(
        self,
        context,
        layout,
        data,
        item,
        icon,
        active_data,
        active_property,
        index=0,
        flt_flag=0,
    ):
        custom_icon = "WORLD"
        if self.layout_type in {"DEFAULT", "COMPACT"}:
            row = layout.row()
            split = row.split(factor=0.1)
            col = split.column()
            col.label(text=f"{index + 1}. ")
            col = split.column()
            col.prop(item, "label", text="", emboss=False)
            hide_icon = "HIDE_OFF" if item.visible else "HIDE_ON"
            row.prop(item, "visible", icon_only=True, icon=hide_icon)
        elif self.layout_type in {"GRID"}:
            layout.alignment = "CENTER"
            layout.label(text="", icon=custom_icon)

    def filter_items(self, context, data, propname):
        if data is None:
            return [], []
        helper_funcs = bpy.types.UI_UL_list
        items = getattr(data, propname)
        # Filter
        filtered = []
        if self.filter_name:
            filtered = helper_funcs.filter_items_by_name(
                self.filter_name,
                self.bitflag_filter_item,
                items,
                "label",
                reverse=self.use_filter_invert,
            )
        if not filtered:
            filtered = [self.bitflag_filter_item] * len(items)
        # Sort
        ordered = []
        if self.use_filter_sort_alpha:
            ordered = helper_funcs.sort_items_by_name(items, "label")
        return filtered, ordered


class MN_PT_Annotations(bpy.types.Panel):
    """
    Panel for annotations
    """

    bl_idname = "MN_PT_annotations"
    bl_label = "Annotations"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Molecular Nodes"

    @classmethod
    def poll(cls, context):
        """Visible only if entity selected is a trajectory or molecule"""
        obj = get_active_entity_object(context)
        if obj is None or context.scene.MNSession.get(obj.uuid) is None:
            return False
        return obj.mn.entity_type in (
            EntityType.MD.value,
            EntityType.MD_STREAMING.value,
            EntityType.MOLECULE.value,
            EntityType.DENSITY.value,
        )

    def draw(self, context):
        scene = context.scene
        object = get_active_entity_object(context)
        if object is None:
            return
        uuid = object.uuid
        entity = scene.MNSession.get(uuid)
        if entity is None:
            return
        annotations_active_index = object.mn.annotations_active_index
        valid_selection = annotations_active_index != -1

        layout = self.layout
        assert layout is not None
        row = layout.row()
        row.prop(object.mn, "annotations_visible", text="Visible")

        row = layout.row()

        MN_UL_AnnotationsList.seqno = 1
        row.template_list(
            "MN_UL_AnnotationsList",
            "annotations_list",
            object,
            "mn_annotations",
            object.mn,
            "annotations_active_index",
            rows=3,
        )
        col = row.column()
        row = col.row()
        op = row.operator("mn.add_annotation", icon="ADD", text="")
        op.uuid = uuid
        row = col.row()
        op = row.operator("mn.remove_annotation", icon="REMOVE", text="")
        if valid_selection:
            op.uuid = uuid
            op.annotation_uuid = object.mn_annotations[annotations_active_index].name
        else:
            row.enabled = False

        if not valid_selection:
            return

        item = object.mn_annotations[annotations_active_index]
        row = layout.row()
        row.prop(item, "type")
        row.enabled = False
        entity_annotation_type = f"{entity._get_annotation_entity_type()}_{item.type}"
        inputs = getattr(item, entity_annotation_type, None)
        instance = entity.annotations._interfaces.get(inputs.uuid)._instance
        if inputs is not None:
            if instance._draw_error is not None:
                row = layout.row()
                row.alert = True
                row.label(text=instance._draw_error, icon="ERROR")

            for prop_name in inputs.__annotations__.keys():
                if prop_name == "uuid":
                    continue
                row = layout.row()
                nbattr = f"_custom_{prop_name}"  # non blender property
                if hasattr(instance, nbattr):
                    # indicate use of non blender property in draw
                    row.label(icon="ERROR")
                    row.alert = True
                if prop_name not in instance._invalid_inputs:
                    row.prop(inputs, prop_name)
                else:
                    row.alert = True
                    row.prop(inputs, prop_name)
                    row = layout.row()
                    row.alert = True
                    row.label(
                        text=instance._invalid_input_messages[prop_name], icon="ERROR"
                    )

        # Add all the common annotation params within the 'Options' panel
        header, panel = layout.panel("annotation_options", default_closed=True)
        header.label(text="Options")

        if panel:
            text_header, text_panel = panel.panel("text_options", default_closed=True)
            text_header.label(text="Text")
            line_header, line_panel = panel.panel("line_options", default_closed=True)
            line_header.label(text="Lines")
            mesh_header, mesh_panel = panel.panel("mesh_options", default_closed=True)
            mesh_header.label(text="Meshes")

            for prop in item.bl_rna.properties:
                if not prop.is_runtime:
                    continue
                if prop.identifier in ("label", "type", "visible"):
                    continue
                if prop.type == "POINTER" and prop.identifier != "mesh_material":
                    continue
                if text_panel and prop.identifier.startswith("text_"):
                    row = text_panel.row()
                    row.prop(item, prop.identifier)
                    if prop.identifier == "text_falloff":
                        row.enabled = item.text_depth
                elif line_panel and prop.identifier.startswith("line_"):
                    row = line_panel.row()
                    row.prop(item, prop.identifier)
                    if prop.identifier == "line_mode":
                        continue
                    row.enabled = item.line_mode != "mesh"
                elif mesh_panel and prop.identifier.startswith("mesh_"):
                    row = mesh_panel.row()
                    row.prop(item, prop.identifier)


class MN_PT_Compositor(bpy.types.Panel):
    """
    Panel for Compositor
    """

    bl_idname = "MN_PT_compositor"
    bl_label = "Options"
    bl_space_type = "NODE_EDITOR"
    bl_region_type = "UI"
    bl_category = "Molecular Nodes"

    @classmethod
    def poll(cls, context):
        return (
            context.area.type == "NODE_EDITOR"
            and context.space_data.tree_type == "CompositorNodeTree"
        )

    def draw(self, context):
        layout = self.layout
        box = layout.box()
        row = box.row()
        row.operator("mn.setup_compositor")


CLASSES = [
    MN_MT_Add,
    MN_MT_Import,
    MN_PT_Scene,
    MN_PT_Object,
    MN_UL_EntitiesList,
    MN_PT_Entities,
    MN_PT_trajectory,
    MN_PT_trajectory_dssp,
    MN_UL_StylesList,
    MN_PT_Styles,
    MN_UL_AnnotationsList,
    MN_PT_Annotations,
    MN_PT_Compositor,
]
