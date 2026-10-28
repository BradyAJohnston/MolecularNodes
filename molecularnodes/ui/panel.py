from typing import cast
import bpy
from bpy.types import UILayout
from ..blender.registered import (
    object_annotations,
    object_props,
    object_selections,
    object_uuid,
    scene_props,
)
from ..entities import Molecule, StreamingTrajectory, molecule
from ..entities.base import EntityType
from ..nodes._utils import get_star_node
from ..session import get_session
from .props import TrajectorySelectionItem
from .utils import check_online_access_for_ui, menu_draw_func, operator_props


def import_options(layout: UILayout) -> None:
    # fetching requires online access, so gate only those entries
    online = check_online_access_for_ui(layout.column())
    op = operator_props(
        online.operator("mn.import_molecule", text="Fetch from PDB", icon="IMPORT")
    )
    op.method = "fetch"
    op.database = "wwpdb"
    op = operator_props(
        online.operator(
            "mn.import_molecule", text="Fetch from AlphaFold", icon="IMPORT"
        )
    )
    op.database = "alphafold"
    op = operator_props(
        online.operator("mn.import_molecule", text="Fetch from PubChem", icon="IMPORT")
    )
    op.method = "fetch"
    op.database = "pubchem"

    op = operator_props(
        layout.operator("mn.import_molecule", text="Import Local File", icon="IMPORT")
    )
    op.method = "local"

    layout.separator()
    op = operator_props(
        layout.operator("mn.import_ensemble", text="Starfile", icon="IMPORT")
    )
    op.ensemble_type = "starfile"
    op = operator_props(
        layout.operator("mn.import_ensemble", text="CellPack", icon="IMPORT")
    )
    op.ensemble_type = "cellpack"

    layout.separator()
    # MD trajectories are loaded through the same unified import dialog
    op = operator_props(
        layout.operator("mn.import_molecule", text="MD Trajectory", icon="IMPORT")
    )
    op.method = "local"
    op = layout.operator("mn.import_oxdna", text="oxDNA", icon="IMPORT")

    layout.separator()
    layout.operator("mn.import_density", text="Density Map", icon="IMPORT")


class MN_MT_Add(bpy.types.Menu):
    bl_idname = "MN_MT_Add"
    bl_label = "Molecular Nodes"

    def draw(self, context) -> None:
        layout = self.layout
        assert layout
        # the Add menu defaults to EXEC_REGION_WIN, which would run the operator
        # directly; force INVOKE so the operators' popup dialogs are shown instead
        with context.temp_override(operator_context="INVOKE_DEFAULT"):
            import_options(layout)


@menu_draw_func
def add_menu_options(self: bpy.types.Menu, context: bpy.types.Context) -> None:
    layout = self.layout
    assert layout
    layout.menu("MN_MT_Add", icon="PARTICLES")


class MN_MT_AddFrameAttribute(bpy.types.Menu):
    """Selections, and calculations and analyses, that can be added to a molecule."""

    bl_idname = "MN_MT_AddFrameAttribute"
    bl_label = "Add"

    def draw(self, context) -> None:
        layout = self.layout
        assert layout
        obj = getattr(context, "object", None)
        entity = get_session(context).get(object_uuid(obj)) if obj is not None else None
        op = operator_props(
            layout.operator(
                "mn.trajectory_selection_add",
                text="Selection",
                icon="RESTRICT_SELECT_OFF",
            )
        )
        if obj is not None:
            op.uuid = object_uuid(obj)
        layout.separator()
        layout.label(text="Calculations")
        row = layout.row()
        # DSSP can be added once, as the `sec_struct` calculation
        row.enabled = (
            isinstance(entity, Molecule) and not entity.calculations._dssp.is_setup
        )
        op = operator_props(row.operator("mn.dssp_init", text="DSSP"))
        if obj is not None:
            op.uuid = object_uuid(obj)
        op.selection = "protein"
        row = layout.row()
        row.enabled = (
            isinstance(entity, Molecule) and not entity.calculations._hbonds.is_setup
        )
        op = operator_props(row.operator("mn.hbonds_init", text="Hydrogen Bonds"))
        if obj is not None:
            op.uuid = object_uuid(obj)
        op.selection = "all"


class MN_MT_Import(bpy.types.Menu):
    bl_idname = "MN_MT_Import"
    bl_label = "Import"

    def draw(self, context) -> None:
        layout = self.layout
        assert layout
        # force INVOKE so the operators' popup dialogs are shown instead of
        # being executed directly
        with context.temp_override(operator_context="INVOKE_DEFAULT"):
            import_options(layout)


@menu_draw_func
def pt_object_context(self, context):
    obj = context.active_object
    if obj is None or not isinstance(get_session(context).match(obj), Molecule):
        return None
    layout = self.layout
    layout.separator()
    layout.operator("mn.align")


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
    obj = context.active_object
    assert obj is not None
    modifier = obj.modifiers["Molecular Nodes"]
    assert isinstance(modifier, bpy.types.NodesModifier)
    ntree = modifier.node_group
    assert ntree is not None
    assert node.node_tree is not None and node.node_tree.interface is not None

    tree = node.node_tree.interface.items_tree

    for item in tree.values():
        if isinstance(item, bpy.types.NodeTreeInterfacePanel):
            col.label(text=item.name)
        elif not isinstance(item, bpy.types.NodeTreeInterfaceSocket):
            continue
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
        op = operator_props(
            box.operator("wm.url_open", text="Selection Langauge Docs", icon="URL")
        )
        op.url = (
            "https://docs.mdanalysis.org/stable/documentation_pages/selections.html"
        )

    return row


def layout_trajectory_playback(
    layout: UILayout, traj: Molecule, panel: bool = True
) -> None:
    if panel:
        header, body = layout.panel(idname="layout_playback")
        header.label(text="Trajectory Playback", icon="OPTIONS")
        if body is None:
            return
        layout = body
    props = traj.props
    is_streaming = isinstance(traj, StreamingTrajectory)

    # a single static structure (one topology, no trajectory) has nothing to play back,
    # so the playback controls are only shown when there are multiple frames
    if not is_streaming and props.n_frames <= 1:
        if panel:
            layout.label(text="Single frame; no playback")
        return

    if is_streaming:
        label = "Streaming trajectory; cannot alter playback"
    else:
        label = "Trajectory has {} {}".format(
            props.n_frames, "frames" if props.n_frames else "frame"
        )

    layout.label(text=label)
    playback = layout.column()
    playback.enabled = not is_streaming
    row = playback.row()

    col = row.column()
    if props.update_with_scene:
        col.prop(props, "frame_hidden")
    else:
        col.prop(props, "frame")
        # with subframes the manual frame is stepped through them, so show the
        # trajectory frame it resolved to
        if props.subframes > 0:
            subrow = col.row()
            subrow.prop(props, "frame_hidden", text="Trajectory Frame")
            subrow.enabled = False
    col.enabled = not props.update_with_scene
    row.prop(props, "update_with_scene")
    row = playback.row()
    col = row.column()
    col.prop(props, "average")
    col.prop(props, "subframes")
    # the offset shifts where playback starts on the timeline, so it only
    # applies when updating with the scene
    subrow = col.row()
    subrow.prop(props, "offset")
    subrow.enabled = props.update_with_scene
    col = row.column()

    # only enable this as an option if the universe has a periodic box
    row = col.row()
    row.prop(props, "correct_periodic")
    row.enabled = traj._has_periodic_box
    col.prop(props, "interpolate")

    # bake frames into a collection for use directly in geometry nodes
    playback.separator()
    playback.operator("mn.frames_to_collection", icon="RENDERLAYERS")


def layout_frame_attributes(
    layout: UILayout, traj: Molecule, panel: bool = True
) -> None:
    """Draw the attributes written on frame changes, with a filter by source, an
    add menu and the selected attribute's details."""
    if panel:
        header, body = layout.panel(idname="layout_frame_attributes")
        header.label(text="Attributes", icon="SPREADSHEET")
        if body is None:
            return
        layout = body
    obj = traj.object
    props = traj.props
    layout.row().prop(props, "frame_attributes_filter", expand=True)
    row = layout.row()
    row.template_list(
        "MN_UL_FrameAttributesList",
        "frame_attributes_list",
        props,
        "frame_attributes",
        props,
        "frame_attributes_index",
        rows=5,
    )
    col = row.column(align=True)
    col.context_pointer_set("object", obj)
    col.menu("MN_MT_AddFrameAttribute", icon="ADD", text="")
    layout_frame_attribute_details(layout, traj)


def layout_frame_attribute_details(layout: UILayout, traj: Molecule) -> None:
    """Draw where the selected attribute comes from, and the settings of a
    selection or calculation."""
    obj = traj.object
    items = traj.props.frame_attributes
    index = traj.props.frame_attributes_index
    if not 0 <= index < len(items):
        return
    item = items[index]
    source_property = cast(bpy.types.EnumProperty, item.bl_rna.properties["source"])
    source = source_property.enum_items[item.source]
    layout.label(text=f"Source: {source.name}")
    if item.source == "SELECTION":
        selection = object_selections(obj).get(item.name)
        if selection is not None:
            selection_string_input(layout, selection)
            row = layout.row()
            row.prop(selection, "updating", text="Update Every Frame")
            row.prop(selection, "periodic")
            row.enabled = not selection.from_atomgroup
        return
    if item.source != "CALCULATION":
        return
    if item.name == "sec_struct" and traj.calculations._dssp.is_setup:
        layout_dssp_settings(layout, traj)
        return
    if item.name == "hbonds" and traj.calculations._hbonds.is_setup:
        layout_hbonds_settings(layout, traj)
        return
    calculation = traj.props.calculations.get(item.name)
    if calculation is not None:
        layout.prop(calculation, "updating", text="Update Every Frame")


def layout_hbonds_settings(layout: UILayout, traj: Molecule) -> None:
    """Draw the settings of the `hbonds` calculation."""
    obj = traj.object
    props = traj.props.hbonds
    layout.prop(
        traj.props.calculations["hbonds"], "updating", text="Update Every Frame"
    )
    row = layout.row()
    row.prop(props, "selection")
    op = operator_props(row.operator("mn.hbonds_init", text="", icon="FILE_REFRESH"))
    op.uuid = object_uuid(obj)
    op.selection = props.selection
    col = layout.column(align=True)
    col.prop(props, "d_a_cutoff")
    col.prop(props, "d_h_a_angle_cutoff")
    col.prop(props, "d_h_cutoff")
    data_object = traj.calculations._hbonds.object
    if data_object is not None:
        mesh = cast(bpy.types.Mesh, data_object.data)
        layout.label(text=f"{len(mesh.edges)} hydrogen bonds")


def layout_dssp_settings(layout: UILayout, traj: Molecule) -> None:
    """Draw the settings of the DSSP `sec_struct` calculation."""
    uuid = object_uuid(traj.object)
    props = traj.props.dssp
    layout.prop(
        traj.props.calculations["sec_struct"], "updating", text="Update Every Frame"
    )
    row = layout.row()
    row.prop(props, "selection")
    op = operator_props(row.operator("mn.dssp_init", text="", icon="FILE_REFRESH"))
    op.uuid = uuid
    op.selection = props.selection
    # display options
    if traj._entity_type == EntityType.MD_STREAMING:
        layout.prop(props, "display_option_streaming")
    else:
        layout.prop(props, "display_option")
    # display option specific params
    if props.display_option == "sliding-window-average":
        layout.prop(props, "window_size")
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
        # averaging the whole trajectory is slow, so it waits for Apply
        row = layout.row()
        split = row.split(factor=0.5)
        col = split.column()
        op = operator_props(col.operator("mn.dssp_apply"))
        op.uuid = uuid
        op.apply_ta_threshold = props.apply_ta_threshold
        op.ta_threshold = props.ta_threshold
        col.enabled = not props.applied
        col = split.column()
        op = operator_props(col.operator("mn.dssp_cancel"))
        op.uuid = uuid


def panel_md_properties(layout, context):
    obj = context.active_object
    session = get_session()
    traj = session.match(obj)
    if traj is not None and not isinstance(traj, molecule.Molecule):
        raise TypeError(f"Expected a trajectory, got {type(traj)}")

    col = layout.column()
    col.enabled = False
    if traj is None:
        col.enabled = True
        col.label(text="Object not linked to a trajectory, please reload one")
        col.prop(obj.mn, "filepath_topology")
        col.prop(obj.mn, "filepath_trajectory")
        col.operator("mn.session_reload_item", text="Reload Trajectory")
        return None

    layout_trajectory_playback(layout, traj)
    layout_frame_attributes(layout, traj)


def panel_object(layout: bpy.types.UILayout, context: bpy.types.Context):
    object = cast(bpy.types.Object, context.active_object)
    if object is None:
        # When an object is deleted, context.ative_object is None
        return
    row = layout.row()
    row.prop(object_props(object), "entity_type")
    row.enabled = False
    try:
        mol_type = object_props(object).entity_type
    except AttributeError:
        return None
    if not object_props(object).is_entity:
        layout.label(text="No MN object selected")
        return None
    if mol_type.startswith("md") or mol_type == EntityType.MOLECULE:
        # molecules and trajectories are both Universe-backed, so both expose selections
        # and (frame-count permitting) playback
        panel_md_properties(layout, context)
    if mol_type in STYLED_ENTITY_TYPES:
        # the same styles list as the viewport sidebar, since styles belong to
        # this object
        header, body = layout.panel(idname="object_styles")
        header.label(text="Styles", icon="MATERIAL")
        if body is not None:
            panel_styles(body, context, object)
    if mol_type == "ensemble-star":
        layout.label(text="Ensemble")
        box = layout.box()
        star_node = get_star_node(object)
        assert isinstance(star_node, bpy.types.GeometryNodeGroup)
        ui_from_node(box, star_node, context=context)
        return None


# must match molecularnodes.scene.world.mn_world_shader_node_name — not imported
# from there because the scene package is an optional dependency of the addon UI
MN_WORLD_SHADER_NODE = "MN_world_shader"


def get_world_shader_node(scene: bpy.types.Scene) -> bpy.types.Node | None:
    """The ``MN_world_shader`` node driving the world background, or None."""
    world = scene.world
    if world is None or world.node_tree is None:
        return None
    return world.node_tree.nodes.get(MN_WORLD_SHADER_NODE)


def panel_render(layout: UILayout, context: bpy.types.Context) -> None:
    """
    The common render settings in one place, mirroring the convenience
    properties that ``mn.Canvas`` exposes for scripting.
    """
    scene = context.scene
    assert scene is not None
    layout.use_property_split = True
    layout.use_property_decorate = False

    col = layout.column()
    col.prop(scene.render, "engine")
    if scene.render.engine == "CYCLES":
        col.prop(scene.cycles, "device")
        col.prop(scene.cycles, "samples")
        col.prop(scene.cycles, "use_denoising", text="Denoise")
    else:
        col.prop(scene.eevee, "taa_render_samples", text="Samples")

    col = layout.column(align=True)
    col.prop(scene.render, "resolution_x", text="Resolution X")
    col.prop(scene.render, "resolution_y", text="Y")
    col.prop(scene.render, "resolution_percentage", text="%")
    layout.prop(scene.render, "film_transparent", text="Transparent")

    col = layout.column()
    # the background shortcuts need the world shader node from the MN preset;
    # sockets match molecularnodes.scene.world.WorldTree
    world_shader = get_world_shader_node(scene)
    if world_shader is not None:
        col.prop(world_shader.inputs[3], "default_value", text="Background")
        col.prop(world_shader.inputs[1], "default_value", text="HDRI Strength")

    header, body = layout.panel(idname="mn_render_camera", default_closed=False)
    header.label(text="Camera")
    if body is not None:
        # aimed at scenes with a single camera, so the scene's active camera
        # is edited directly rather than offering a picker
        camera = scene.camera
        if camera is None:
            body.label(text="No camera in the scene")
        else:
            cam_data = cast(bpy.types.Camera, camera.data)
            dof = cam_data.dof
            assert dof is not None
            body.prop(cam_data, "lens", text="Focal Length")
            body.prop(dof, "use_dof", text="Depth of Field")
            focus = body.column()
            focus.enabled = dof.use_dof
            focus.prop(dof, "focus_object")
            distance = focus.row()
            distance.active = dof.focus_object is None
            row = distance.row(align=True)
            row.prop(dof, "focus_distance", text="Focus Distance")
            # the depth eyedropper samples through the camera, so its poll only
            # passes while the viewport is in camera view — grey it out
            # deliberately (rather than letting the poll flicker it on hover)
            # and hint at how to enable it
            region_3d = getattr(context.space_data, "region_3d", None)
            in_camera_view = (
                region_3d is not None and region_3d.view_perspective == "CAMERA"
            )
            pick = row.row(align=True)
            pick.enabled = in_camera_view
            op = operator_props(
                pick.operator("ui.eyedropper_depth", icon="EYEDROPPER", text="")
            )
            op.prop_data_path = "scene.camera.data.dof.focus_distance"
            if not in_camera_view:
                focus.label(
                    text="Enter camera view to pick focus distance", icon="INFO"
                )
            focus.prop(dof, "aperture_fstop", text="F-Stop")

    header, body = layout.panel(idname="mn_render_color", default_closed=True)
    header.label(text="Color Management")
    if body is not None:
        view = scene.view_settings
        body.prop(view, "view_transform")
        body.prop(view, "look")
        body.prop(view, "exposure")
        body.prop(view, "gamma")

    header, body = layout.panel(idname="mn_render_animation", default_closed=True)
    header.label(text="Animation")
    if body is not None:
        col = body.column(align=True)
        col.prop(scene, "frame_start", text="Frame Start")
        col.prop(scene, "frame_end", text="End")
        body.prop(scene.render, "fps", text="FPS")


class MN_PT_Render(bpy.types.Panel):
    """
    Panel collecting the common render settings (Viewport)
    """

    bl_idname = "MN_PT_render"
    bl_label = "Render"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Molecular Nodes"

    def draw(self, context):
        layout = self.layout
        assert layout
        panel_render(layout, context)


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
    if context.scene is None:
        return None
    index = scene_props(context.scene).entities_active_index
    objects = bpy.data.objects
    if 0 <= index < len(objects):
        obj = objects[index]
        if object_props(obj).is_entity:
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
            linked = get_session(context).get(object_uuid(obj)) is not None
            row.label(text="", icon="LINKED" if linked else "UNLINKED")
            row.prop(obj, "name", text="", emboss=False)
            # use the object viewport visibility to determine the icon
            # we do not have direct callbacks for raw object visibility changes
            mn = object_props(obj)
            hide_icon = "HIDE_OFF" if mn.visible else "HIDE_ON"
            row.prop(mn, "visible", icon_only=True, icon=hide_icon)
        elif self.layout_type in {"GRID"}:
            layout.alignment = "CENTER"
            layout.label(text="", icon="OBJECT_DATA")

    def filter_items(self, context, data, property):
        if data is None:
            return [], []
        objects = getattr(data, property)
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


# the node-editor socket dot for each attribute data type
ATTRIBUTE_TYPE_ICONS = {
    "FLOAT": "NODE_SOCKET_FLOAT",
    "INT": "NODE_SOCKET_INT",
    "INT8": "NODE_SOCKET_INT",
    "INT16_2D": "NODE_SOCKET_INT_VECTOR",
    "INT32_2D": "NODE_SOCKET_INT_VECTOR",
    "BOOLEAN": "NODE_SOCKET_BOOLEAN",
    "FLOAT_VECTOR": "NODE_SOCKET_VECTOR",
    "FLOAT2": "NODE_SOCKET_VECTOR",
    "FLOAT4": "NODE_SOCKET_VECTOR",
    "FLOAT_COLOR": "NODE_SOCKET_RGBA",
    "BYTE_COLOR": "NODE_SOCKET_RGBA",
    "QUATERNION": "NODE_SOCKET_ROTATION",
    "FLOAT4X4": "NODE_SOCKET_MATRIX",
    "STRING": "NODE_SOCKET_STRING",
}


class MN_UL_FrameAttributesList(bpy.types.UIList):
    """
    UIList of the attributes written on frame changes.

    Lists the items of `obj.mn.frame_attributes`, showing each attribute's name
    and the data type of the mesh attribute.
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
        # `data` is the object's `mn` properties
        obj = cast(bpy.types.Object, data.id_data)
        attribute = cast(bpy.types.Mesh, obj.data).attributes.get(item.name)
        if self.layout_type in {"DEFAULT", "COMPACT"}:
            row = layout.row()
            source = item.bl_rna.properties["source"].enum_items[item.source]
            row.label(text=item.name, icon=source.icon)
            calculation = (
                object_props(obj).calculations.get(item.name)
                if item.source == "CALCULATION"
                else None
            )
            if calculation is not None and calculation.object is not None:
                # a calculation that writes an object rather than an attribute
                row.label(text="Object", icon="OBJECT_DATA")
            elif attribute is None:
                # not stored yet, or a selection that failed
                row.label(text="")
            else:
                row.label(
                    text=attribute.bl_rna.properties["data_type"]
                    .enum_items[attribute.data_type]
                    .name,
                    icon=ATTRIBUTE_TYPE_ICONS.get(attribute.data_type, "NONE"),
                )
            # selections and calculations can be set to not update per frame
            if item.source == "SELECTION":
                settings = object_selections(obj).get(item.name)
            elif item.source == "CALCULATION":
                settings = object_props(obj).calculations.get(item.name)
            else:
                settings = None
            if settings is not None:
                row.prop(settings, "updating", text="", icon="FILE_REFRESH")
            if item.source != "TRAJECTORY":
                op = operator_props(
                    row.operator(
                        "mn.frame_attribute_remove", text="", icon="X", emboss=False
                    )
                )
                op.uuid = object_uuid(obj)
                op.name = item.name
        elif self.layout_type in {"GRID"}:
            layout.alignment = "CENTER"
            layout.label(text=item.name)

    def filter_items(self, context, data, property):
        """Show the attributes from the source picked above the list, and apply
        the list's own name filter and sorting."""
        items = getattr(data, property)
        helper = bpy.types.UI_UL_list
        if self.filter_name:
            # the stubs type the helper as returning None, rather than the flags
            flags = cast(
                list[int],
                helper.filter_items_by_name(
                    self.filter_name, self.bitflag_filter_item, items, "name"
                ),
            )
        else:
            flags = [self.bitflag_filter_item] * len(items)
        source = data.frame_attributes_filter
        if source != "ALL":
            flags = [
                flag if item.source == source else 0 for flag, item in zip(flags, items)
            ]
        order = (
            helper.sort_items_by_name(items, "name")
            if self.use_filter_sort_alpha
            else []
        )
        return flags, order


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
        entity = get_session(context).get(object_uuid(obj)) if obj is not None else None

        col = row.column()
        # reload/relink the entity into the session from its recorded source
        row = col.row()
        row.operator("mn.session_reload_item", icon="FILE_REFRESH", text="")
        row = col.row()
        op = operator_props(
            row.operator("mn.session_remove_item", icon="REMOVE", text="")
        )
        if entity is None or obj is None:
            row.enabled = False
        else:
            op.uuid = object_uuid(obj)

        if obj is None:
            return
        # display details of the selected entity object
        row = layout.row()
        row.prop(object_props(obj), "entity_type")
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
        if obj is None or get_session(context).get(object_uuid(obj)) is None:
            return False
        mn = object_props(obj)
        if mn.entity_type not in (
            EntityType.MD_STREAMING,
            EntityType.MD_OXDNA,
            EntityType.MOLECULE,
        ):
            return False
        # a single static structure has no playback; streaming has an unknown frame count
        return mn.entity_type == EntityType.MD_STREAMING or mn.n_frames > 1

    def draw(self, context):
        layout = cast(UILayout, self.layout)
        assert layout
        # To enable the animatate dot next to property in UI
        # layout.use_property_split = True
        # layout.use_property_decorate = True
        obj = get_active_entity_object(context)
        # poll ensures a linked trajectory is active
        assert obj is not None
        traj = cast(Molecule, get_session(context).get(object_uuid(obj)))

        layout_trajectory_playback(layout, traj, panel=False)


class MN_PT_Attributes(bpy.types.Panel):
    """
    Panel for the attributes a molecule writes on frame changes: from the
    trajectory, its selections and its calculations
    """

    bl_idname = "MN_PT_attributes"
    bl_label = "Attributes"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Molecular Nodes"

    @classmethod
    def poll(cls, context):
        """Visible for a Universe-backed entity, with or without playback."""
        obj = get_active_entity_object(context)
        if obj is None or obj.type != "MESH":
            return False
        if get_session(context).get(object_uuid(obj)) is None:
            return False
        return object_props(obj).entity_type in (
            EntityType.MD_STREAMING,
            EntityType.MD_OXDNA,
            EntityType.MOLECULE,
        )

    def draw(self, context):
        layout = cast(UILayout, self.layout)
        obj = get_active_entity_object(context)
        # poll ensures a linked trajectory is active
        assert obj is not None
        traj = cast(Molecule, get_session(context).get(object_uuid(obj)))
        layout_frame_attributes(layout, traj, panel=False)


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
                input = cast(bpy.types.NodeSocketBool, item.inputs["Visible"])
                hide_icon = "HIDE_OFF" if input.default_value else "HIDE_ON"
                row.prop(input, "default_value", icon_only=True, icon=hide_icon)
        elif self.layout_type in {"GRID"}:
            layout.alignment = "CENTER"
            layout.label(text="", icon="WORLD")

    def filter_items(self, context, data, property):
        if data is None:
            return [], []
        items = getattr(data, property)
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


# entity types whose objects carry styles in a "Molecular Nodes" modifier tree
STYLED_ENTITY_TYPES = (
    EntityType.MD_STREAMING,
    EntityType.MOLECULE,
    EntityType.DENSITY,
)


def panel_styles(
    layout: UILayout, context: bpy.types.Context, obj: bpy.types.Object
) -> None:
    """
    The styles list and controls for an object, shared between the viewport
    Styles panel and the Object properties panel.
    """
    # style nodes live in the object's node tree — read them directly from
    # Blender data so the panel does not depend on the session being linked
    node_group = get_entity_node_group(obj)
    if node_group is None:
        layout.label(text="No Molecular Nodes modifier on this object")
        return

    # the style node selected in the list, if any
    mn = object_props(obj)
    index = mn.styles_active_index
    style_node = None
    if 0 <= index < len(node_group.nodes):
        node = node_group.nodes[index]
        if is_style_node(node):
            style_node = cast(bpy.types.GeometryNodeGroup, node)

    entity = get_session(context).get(object_uuid(obj))

    # list the style nodes in the tree and let the user select one
    row = layout.row()
    row.template_list(
        "MN_UL_StylesList",
        "styles_list",
        node_group,
        "nodes",
        mn,
        "styles_active_index",
        rows=3,
    )
    col = row.column()
    # adding uses the entity API when linked, and otherwise builds the
    # style branch directly in the object's node tree; density styles are
    # different nodes, so density objects only get swap/remove
    add = col.row()
    add.enabled = isinstance(entity, Molecule) or (
        entity is None and mn.entity_type != EntityType.DENSITY
    )
    op = operator_props(add.operator("mn.add_style", icon="ADD", text=""))
    op.uuid = object_uuid(obj)
    op.name_object = obj.name
    remove = col.row()
    remove.enabled = style_node is not None
    op = operator_props(remove.operator("mn.remove_style", icon="REMOVE", text=""))
    if style_node is not None:
        op.name_tree = node_group.name
        op.name_node = style_node.name

    if style_node is None:
        layout.label(text="Select a style to edit its properties")
        return

    # swap the selected style node for a different style
    assert style_node.node_tree is not None
    row = layout.row(align=True)
    row.label(text="Style:")
    op = operator_props(
        row.operator_menu_enum(
            "mn.swap_style",
            "style",
            text=style_node.node_tree.name.replace("Style ", ""),
        )
    )
    op.name_tree = node_group.name
    op.name_node = style_node.name

    # display the selection string if using a named attribute
    if style_node.inputs["Selection"].links and entity is not None:
        node = style_node.inputs["Selection"].links[0].from_node
        if isinstance(node, bpy.types.GeometryNodeInputNamedAttribute):
            if isinstance(entity, Molecule):
                assert node.inputs is not None
                selection = entity.selections.get(node.inputs["Name"].default_value)
                layout.prop(selection, "string", text="Selection")
    else:
        op = operator_props(layout.operator(operator="mn.add_selection_to_style"))
        op.node_tree = node_group.name
        op.node_name = style_node.name

    # display the selected style node's name and its input properties
    header, body = layout.panel(idname="style_properties")
    header.label(text=style_node.label or style_node.name)
    if body is not None:
        body.template_node_inputs(style_node)


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
        return obj is not None and object_props(obj).entity_type in STYLED_ENTITY_TYPES

    def draw(self, context):
        layout = self.layout
        assert layout is not None
        obj = get_active_entity_object(context)
        if obj is None:
            return
        panel_styles(layout, context, obj)


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

    def filter_items(self, context, data, property):
        if data is None:
            return [], []
        helper_funcs = bpy.types.UI_UL_list
        items = getattr(data, property)
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
    bl_options = {"DEFAULT_CLOSED"}

    @classmethod
    def poll(cls, context):
        """Visible only if entity selected is a trajectory or molecule"""
        obj = get_active_entity_object(context)
        if obj is None or get_session(context).get(object_uuid(obj)) is None:
            return False
        return object_props(obj).entity_type in (
            EntityType.MD_STREAMING,
            EntityType.MOLECULE,
            EntityType.DENSITY,
        )

    def draw(self, context):
        scene = context.scene
        object = get_active_entity_object(context)
        if object is None:
            return
        uuid = object_uuid(object)
        entity = scene.MNSession.get(uuid)
        if entity is None:
            return
        mn = object_props(object)
        annotations = object_annotations(object)
        annotations_active_index = mn.annotations_active_index
        valid_selection = annotations_active_index != -1

        layout = self.layout
        assert layout is not None
        row = layout.row()
        row.prop(mn, "annotations_visible", text="Visible")

        row = layout.row()

        row.template_list(
            "MN_UL_AnnotationsList",
            "annotations_list",
            object,
            "mn_annotations",
            mn,
            "annotations_active_index",
            rows=3,
        )
        col = row.column()
        row = col.row()
        op = operator_props(row.operator("mn.add_annotation", icon="ADD", text=""))
        op.uuid = uuid
        row = col.row()
        op = operator_props(
            row.operator("mn.remove_annotation", icon="REMOVE", text="")
        )
        if valid_selection:
            op.uuid = uuid
            op.annotation_uuid = annotations[annotations_active_index].name
        else:
            row.enabled = False

        if not valid_selection:
            return

        item = annotations[annotations_active_index]
        row = layout.row()
        row.prop(item, "type")
        row.enabled = False
        entity_annotation_type = f"{entity._get_annotation_entity_type()}_{item.type}"
        inputs = getattr(item, entity_annotation_type, None)
        if inputs is not None:
            instance = entity.annotations._interfaces.get(inputs.uuid)._instance
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
        assert layout is not None
        box = layout.box()
        row = box.row()
        row.operator("mn.setup_compositor")


CLASSES = [
    MN_MT_Add,
    MN_MT_AddFrameAttribute,
    MN_MT_Import,
    MN_PT_Scene,
    MN_PT_Object,
    MN_UL_EntitiesList,
    MN_PT_Entities,
    MN_UL_FrameAttributesList,
    MN_PT_trajectory,
    MN_PT_Attributes,
    MN_UL_StylesList,
    MN_PT_Styles,
    MN_UL_AnnotationsList,
    MN_PT_Annotations,
    MN_PT_Render,
    MN_PT_Compositor,
]
