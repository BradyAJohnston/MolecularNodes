import bpy
from bpy.props import (
    BoolProperty,
    CollectionProperty,
    EnumProperty,
    FloatProperty,
    IntProperty,
    PointerProperty,
    StringProperty,
)
from databpy.object import LinkedObjectError
from ..blender.utils import set_object_visibility
from ..entities.base import EntityType
from ..handlers import _update_entities
from ..session import get_entity

uuid_property = StringProperty(
    name="UUID",
    description="Unique ID for referencing the required objects in the MNSession",
    default="",
)

ENTITY_ITEMS = (
    ("None", "None", "Not an MN entity"),
    ("molecule", "Molecule", "A single molecule"),
    ("ensemble", "Ensemble", "A collection of molecules"),
    ("density", "Density", "A density grid"),
    ("md", "Trajectory", "A molecular dynamics trajectory"),
    ("md-oxdna", "oxDNA Trajectory", "A oxDNA molecular dynamics trajectory"),
    (
        "md-streaming",
        "Streaming Trajectory",
        "A streaming IMD molecular dynamics trajectory",
    ),
    ("ensemble-star", "Star Ensemble", "A starfile ensemble"),
    ("ensemble-cellpack", "CellPack Ensemble", "A CellPack model ensemble"),
)

SURFACE_STYLE_ITEMS = (
    (
        "density_surface",
        "Surface",
        "Style Density Surface",
    ),
    # (
    #     "density_iso_surface",
    #     "Iso Surface",
    #     "Style Density ISO Surface",
    # ),
    (
        "density_wire",
        "Wire",
        "Style Density Wire",
    ),
)


def _get_frame(self):
    return self.get("frame", 0)


def _set_frame(self, frame):
    if frame >= self.n_frames:
        frame = self.n_frames - 1
    self["frame"] = frame
    _update_entities(self, bpy.context)


def _get_entity_visibility(self) -> bool:
    """get callback for entity visibility property"""
    return self.get("visible", True)


def _set_entity_visibility(self, visible: bool) -> None:
    """set callback for entity visibility property"""
    self["visible"] = visible
    entity = bpy.context.scene.MNSession.get(self.name)
    if entity is not None:
        set_object_visibility(entity.object, self.visible)
        entity.annotations._update_annotation_object()


def _entities_active_index_callback(self, context: bpy.context) -> None:  # type: ignore
    """update callback for entities active_index change"""
    if self.entities_active_index == -1:
        return
    uuid = context.scene.mn.entities[self.entities_active_index].name
    try:
        # object might not yet be created during session entity registration
        entity_object = context.scene.MNSession.get(uuid).object
    except (LinkedObjectError, AttributeError):
        return
    # just setting view_layer.objects.active is not enough
    bpy.ops.object.select_all(action="DESELECT")  # deselect all objects
    if entity_object.name in context.view_layer.objects:
        context.view_layer.objects.active = entity_object  # make active object
    bpy.context.view_layer.update()  # update view layer to reflect changes
    if bpy.context.active_object:  # can be None for hidden objects
        bpy.context.active_object.select_set(True)  # set as selected object


class EntityProperties(bpy.types.PropertyGroup):
    # name property is implicit and is set to uuid for find lookups
    # type value is one of EntityType enum
    __slots__ = []
    type: StringProperty(name="Entity Type", default="")  # type: ignore
    visible: BoolProperty(  # type: ignore
        name="visible",
        description="Visibility of the entity",
        default=True,
        get=_get_entity_visibility,
        set=_set_entity_visibility,
    )  # type: ignore


def _update_dssp_display_option(self, context):
    entity = context.scene.MNSession.get(self.id_data.uuid)
    if entity is None or self.cancelling:
        return
    if entity._entity_type == EntityType.MD_STREAMING:
        display_option = getattr(self, "display_option_streaming")
    else:
        display_option = getattr(self, "display_option")
    # call none and per-frame directly
    if display_option == "none":
        entity.dssp.show_none()
        _update_entities(self, context)
    elif display_option == "per-frame":
        entity.dssp.show_per_frame()
        _update_entities(self, context)
    elif display_option == "sliding-window-average":
        sw_threshold = self.sw_threshold if self.apply_sw_threshold else None
        entity.dssp.show_sliding_window_average(self.window_size, sw_threshold)
        _update_entities(self, context)
    else:
        self.applied = False


def _update_dssp_applied(self, context):
    if self.applied:
        _update_entities(self, context)


class DSSPProperties(bpy.types.PropertyGroup):
    display_option: EnumProperty(  # type: ignore
        name="Display",
        description="Options to display secondary structures",
        items=(
            ("none", "None", "Do not show secondary structures"),
            ("per-frame", "Per Frame", "Secondary structures calculated per frame"),
            (
                "sliding-window-average",
                "Sliding Window Average",
                "Average secondary structures of a sliding window of frames",
            ),
            (
                "trajectory-average",
                "Trajectory Average",
                "Average secondary structures across all frames",
            ),
        ),
        default="per-frame",
        update=_update_dssp_display_option,
    )
    display_option_streaming: EnumProperty(  # type: ignore
        name="Display",
        description="Options to display secondary structures",
        items=(
            ("none", "None", "Do not show secondary structures"),
            ("per-frame", "Per Frame", "Secondary structures calculated per frame"),
        ),
        default="per-frame",
        update=_update_dssp_display_option,
    )
    window_size: IntProperty(  # type: ignore
        name="Window Size",
        description="Number of frames in the sliding window",
        min=1,
        soft_max=10,
        default=5,
        update=_update_dssp_display_option,
    )  # type: ignore
    apply_sw_threshold: BoolProperty(  # type: ignore
        name="Apply Threshold",
        description="Apply a threshold comparison to calculated mean",
        default=False,
        update=_update_dssp_display_option,
    )  # type: ignore
    sw_threshold: FloatProperty(  # type: ignore
        name="Threshold",
        description="Threshold fraction of frames for sliding window average",
        subtype="FACTOR",
        min=0.0,
        max=1.0,
        default=0.5,
        update=_update_dssp_display_option,
    )  # type: ignore
    apply_ta_threshold: BoolProperty(  # type: ignore
        name="Apply Threshold",
        description="Apply a threshold comparison to calculated mean",
        default=False,
        update=_update_dssp_display_option,
    )  # type: ignore
    ta_threshold: FloatProperty(  # type: ignore
        name="Threshold",
        description="Threshold fraction of frames for trajectory average",
        subtype="FACTOR",
        min=0.0,
        max=1.0,
        default=0.5,
        update=_update_dssp_display_option,
    )
    applied: BoolProperty(  # type: ignore
        default=True,
        update=_update_dssp_applied,
    )
    cancelling: BoolProperty(default=False)  # type: ignore


class MolecularNodesSceneProperties(bpy.types.PropertyGroup):
    __slots__ = []
    entities: CollectionProperty(name="Entities", type=EntityProperties)  # type: ignore
    entities_active_index: IntProperty(  # type: ignore
        name="Active entity index",
        default=-1,
        update=_entities_active_index_callback,
    )  # type: ignore

    is_updating: BoolProperty(  # type: ignore
        name="Updating",
        description="Currently updating data in the scene, don't trigger more updates",
        default=False,
    )


def _update_annotations_visibility(self, context):
    entity = context.scene.MNSession.get(self.id_data.uuid)
    if entity is not None:
        if self.annotations_visible:
            entity.annotations._draw_handler_add()
        else:
            entity.annotations._draw_handler_remove()
        entity.annotations._update_annotation_object()


class MolecularNodesObjectProperties(bpy.types.PropertyGroup):
    __slots__ = []
    styles_active_index: IntProperty(default=-1)  # type: ignore
    annotations_active_index: IntProperty(default=-1)  # type: ignore
    annotations_next_index: IntProperty(default=0)  # type: ignore

    annotations_visible: BoolProperty(  # type: ignore
        name="Visible",
        description="Visibility of all annotations",
        default=True,
        update=_update_annotations_visibility,
    )

    biological_assemblies: StringProperty(  # type: ignore
        name="Biological Assemblies",
        description="A list of biological assemblies to be created",
        default="",
    )

    entity_type: EnumProperty(  # type: ignore
        name="Entity Type",
        description="How the file was imported, dictating how MN interacts with it",
        items=ENTITY_ITEMS,
        default="None",
    )

    code: StringProperty(  # type: ignore
        name="PDB",
        description="PDB code used to download this structure",
        maxlen=4,
        options={"HIDDEN"},
    )
    trajectory_selection_index: IntProperty(  # type: ignore
        name="Index of selection",
        description="Index of selection, that is selected for the UI",
        default=0,
    )
    frame_hidden: IntProperty(  # type: ignore
        name="Frame",
        description="Frame of the loaded trajectory",
        default=0,
        min=0,
    )
    frame: IntProperty(  # type: ignore
        name="Frame",
        description="Frame of the loaded trajectory",
        default=0,
        min=0,
        set=_set_frame,
        get=_get_frame,
    )
    n_frames: IntProperty(  # type: ignore
        name="Number of Frames",
        description="Number of frames in the loaded trajectory",
        default=0,
        min=0,
    )
    update_with_scene: BoolProperty(  # type: ignore
        name="Update with Scene",
        description="Update the trajectory with the scene frame",
        default=True,
        update=_update_entities,
    )
    subframes: IntProperty(  # type: ignore
        name="Subframes",
        description="Number of subframes to insert between frames of the loaded trajectory",
        default=0,
        update=_update_entities,
        min=0,
    )
    offset: IntProperty(  # type: ignore
        name="Offset",
        description="Offset the starting playback for the trajectory on the timeine. Positive starts the playback later than frame 0, negative starts it earlier than frame 0",
        default=0,
        update=_update_entities,
    )

    interpolate: BoolProperty(  # type: ignore
        name="Interpolate",
        description="Whether to interpolate when using subframes",
        default=True,
        update=_update_entities,
    )
    average: IntProperty(  # type: ignore
        name="Average",
        description="Average the position this number of frames either side of the current frame",
        default=0,
        update=_update_entities,
        min=0,
        soft_max=5,
    )
    correct_periodic: BoolProperty(  # type: ignore
        name="Correct",
        description="Correct for periodic boundary crossing when using interpolation or averaging. Assumes cubic dimensions and only works if the unit cell is orthorhombic",
        default=False,
        update=_update_entities,
    )
    filepath_trajectory: StringProperty(  # type: ignore
        name="Trajectory",
        description="Filepath for the `trajectory` of the Object",
        subtype="FILE_PATH",
        default="",
    )
    filepath_topology: StringProperty(  # type: ignore
        name="Topology",
        description="Filepath for the Topology of the Object",
        subtype="FILE_PATH",
        default="",
    )
    dssp: PointerProperty(type=DSSPProperties)  # type: ignore


class TrajectorySelectionItem(bpy.types.PropertyGroup):
    """Group of properties for custom selections for MDAnalysis import."""

    __slots__ = []

    name: StringProperty(  # type: ignore
        name="Attribute Name",
        description="Name of the attribute that will be created when storing on the mesh",
    )

    string: StringProperty(  # type: ignore
        name="Selection",
        description="Selection to be applied, written in the MDAnalysis selection language",
        default="name CA",
        update=_update_entities,
    )

    previous_string: StringProperty()  # type: ignore

    updating: BoolProperty(  # type: ignore
        name="Updating",
        description="Potential recalculate the selection when the scene frame changes",
        default=True,
        update=_update_entities,
    )

    periodic: BoolProperty(  # type: ignore
        name="Periodic",
        description="For geometric selections, whether to account for atoms in different periodic images when searching",
        default=True,
        update=_update_entities,
    )

    message: StringProperty(  # type: ignore
        name="Message",
        description="Message to report back from `universe.select_atoms()`",
        default="",
    )

    from_atomgroup: BoolProperty(  # type: ignore
        name="From AtomGroup",
        description="If the UI item has been created from an existing AtomGroup. Will prevent editing in the UI by a user.",
        default=False,
    )


class MN_UL_TrajectorySelectionListUI(bpy.types.UIList):
    """UI List"""

    def draw_item(
        self,
        context,
        layout,
        data,
        item,
        icon,
        active_data,
        active_property,
        *,
        index=0,
        flt_flag=0,
    ):
        custom_icon = "VIS_SEL_11"

        if self.layout_type in {"DEFAULT", "COMPACT"}:
            row = layout.row()
            if item.message != "":
                custom_icon = "ERROR"
                row.alert = True

            col = row.column()
            col.prop(item, "name", text="", emboss=False)
            col.enabled = False
            row.prop(item, "updating", icon_only=True, icon="FILE_REFRESH")
            row.prop(item, "periodic", icon_only=True, icon="CUBE")
            if item.from_atomgroup:
                row.enabled = False

        elif self.layout_type in {"GRID"}:
            layout.alignment = "CENTER"
            layout.label(text="", icon=custom_icon)


class MN_OT_Universe_Selection_Add(bpy.types.Operator):
    "Add a new custom selection to a trajectory"

    bl_idname = "mn.trajectory_selection_add"
    bl_label = "+"
    bl_description = "Add a new boolean attribute for the given MDA selection string"

    def execute(self, context):
        traj = get_entity(context)
        traj.selections.from_string("all")
        traj.selections.ui_index = max(0, len(traj.selections) - 1)
        return {"FINISHED"}


class MN_OT_Universe_Selection_Delete(bpy.types.Operator):
    bl_idname = "mn.trajectory_selection_remove"
    bl_label = "-"
    bl_description = "Delete the given boolean selection from the universe"

    @classmethod
    def poll(cls, context):
        return context.active_object.mn_trajectory_selections

    def execute(self, context):
        traj = get_entity(context)
        index = traj.selections.ui_index
        traj.selections.remove(index)

        # the length of items in the list has changed, set the currently selected index
        # to a new value. Ensure it is between 0 and the length of the items in the list
        traj.selections.ui_index = max(0, min(index, len(traj.selections.ui_items) - 1))

        return {"FINISHED"}


CLASSES = [
    EntityProperties,
    DSSPProperties,
    MolecularNodesObjectProperties,
    MolecularNodesSceneProperties,
    TrajectorySelectionItem,  # item has to be registered the ListUI and to work properly
    MN_UL_TrajectorySelectionListUI,
    MN_OT_Universe_Selection_Add,
    MN_OT_Universe_Selection_Delete,
]
