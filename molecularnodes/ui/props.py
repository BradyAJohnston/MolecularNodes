import json
import bpy
from bpy.props import (
    BoolProperty,
    EnumProperty,
    FloatProperty,
    IntProperty,
    PointerProperty,
    StringProperty,
)
from nodebpy.nodes.geometry import NamedAttribute
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
        "Density Style Surface",
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


def _get_object_visibility(self) -> bool:
    """get callback for object visibility property"""
    try:
        return self.id_data.visible_get()
    except RuntimeError:
        # object not in the current view layer
        return True


def _set_object_visibility(self, visible: bool) -> None:
    """set callback for object visibility property"""
    obj = self.id_data
    set_object_visibility(obj, visible)
    entity = bpy.context.scene.MNSession.get(obj.uuid)
    if entity is not None:
        entity.annotations._update_annotation_object()


def _get_entities_active_index(self) -> int:
    """
    Derive the list selection from the scene's active object.

    The active object is the single source of truth for the Entities list, so
    selecting an entity object in the outliner or viewport is mirrored in the
    list. Returns -1 when the active object is not a molecular entity.
    """
    obj = bpy.context.view_layer.objects.active
    if obj is None or not obj.mn.is_entity:
        return -1
    return bpy.data.objects.find(obj.name)


def _set_entities_active_index(self, index: int) -> None:
    """set callback: mirror the list selection into the scene selection"""
    # the index points into bpy.data.objects (the Entities list data source)
    if index < 0 or index >= len(bpy.data.objects):
        return
    entity_object = bpy.data.objects[index]
    if not entity_object.mn.is_entity:
        return
    context = bpy.context
    if entity_object.name not in context.view_layer.objects:
        return
    for obj in list(context.selected_objects):
        obj.select_set(False)
    context.view_layer.objects.active = entity_object
    entity_object.select_set(True)


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
    )
    apply_sw_threshold: BoolProperty(  # type: ignore
        name="Apply Threshold",
        description="Apply a threshold comparison to calculated mean",
        default=False,
        update=_update_dssp_display_option,
    )
    sw_threshold: FloatProperty(  # type: ignore
        name="Threshold",
        description="Threshold fraction of frames for sliding window average",
        subtype="FACTOR",
        min=0.0,
        max=1.0,
        default=0.5,
        update=_update_dssp_display_option,
    )
    apply_ta_threshold: BoolProperty(  # type: ignore
        name="Apply Threshold",
        description="Apply a threshold comparison to calculated mean",
        default=False,
        update=_update_dssp_display_option,
    )
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
    entities_active_index: IntProperty(  # type: ignore
        name="Active entity index",
        description="Index into bpy.data.objects of the entity object active in"
        " the Entities list, derived from the scene's active object",
        default=-1,
        get=_get_entities_active_index,
        set=_set_entities_active_index,
    )

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

    internal_biological_assemblies: StringProperty(  # type: ignore
        name="Biological Assemblies",
        description="A list of biological assemblies to be created",
        default="",
        options={"HIDDEN"},
    )

    @property
    def biological_assemblies(self) -> dict:
        """Return the biological assemblies for the entity object"""
        if not self.internal_biological_assemblies:
            return {}
        return json.loads(self.internal_biological_assemblies)

    @biological_assemblies.setter
    def biological_assemblies(self, value: dict | str | None) -> None:
        """Set the biological assemblies for the entity object"""
        if value is None or value == "":
            self.internal_biological_assemblies = ""
        elif isinstance(value, str):
            self.internal_biological_assemblies = value
        else:
            self.internal_biological_assemblies = json.dumps(value)

    entity_type: EnumProperty(  # type: ignore
        name="Entity Type",
        description="How the file was imported, dictating how MN interacts with it",
        items=ENTITY_ITEMS,
        default="None",
    )

    visible: BoolProperty(  # type: ignore
        name="Visible",
        description="Visibility of the entity object in the viewport and renders",
        get=_get_object_visibility,
        set=_set_object_visibility,
    )

    code: StringProperty(  # type: ignore
        name="PDB",
        description="PDB code used to download this structure",
        maxlen=12,
        options={"HIDDEN"},
    )
    filepath: StringProperty(  # type: ignore
        name="File Path",
        description="Source file this entity was loaded from, used to reload it",
        subtype="FILE_PATH",
        default="",
        options={"HIDDEN"},
    )
    database: StringProperty(  # type: ignore
        name="Database",
        description="Database this structure was fetched from, used to reload it",
        default="",
        options={"HIDDEN"},
    )

    internal_chain_ids: StringProperty(  # type: ignore
        name="Chain IDs",
        description="A list of chain IDs for the entity object",
        default="",
        options={"HIDDEN"},
    )

    @property
    def chain_ids(self) -> list[str]:
        """Return a list of chain IDs for the entity object"""
        return self.internal_chain_ids.split(",") if self.internal_chain_ids else []

    @chain_ids.setter
    def chain_ids(self, value: list[str] | None) -> None:
        """Set the list of chain IDs for the entity object"""
        if value is None:
            self.internal_chain_ids = ""
        else:
            self.internal_chain_ids = ",".join(value)

    internal_entity_ids: StringProperty(  # type: ignore
        name="Entity IDs",
        description="A list of entity IDs for the entity object",
        default="",
        options={"HIDDEN"},
    )

    @property
    def entity_ids(self) -> list[str]:
        """Return a list of entity IDs for the entity object"""
        return self.internal_entity_ids.split(",") if self.internal_entity_ids else []

    @entity_ids.setter
    def entity_ids(self, value: list[str] | None) -> None:
        """Set the list of entity IDs for the entity object"""
        if value is None:
            self.internal_entity_ids = ""
        else:
            self.internal_entity_ids = ",".join(value)

    internal_segments: StringProperty(  # type: ignore
        name="Segments",
        description="A list of segment IDs for the entity object",
        default="",
        options={"HIDDEN"},
    )

    @property
    def segments(self) -> list[str]:
        """Return a list of segment IDs for the entity object"""
        return self.internal_segments.split(",") if self.internal_segments else []

    @segments.setter
    def segments(self, value: list[str] | None) -> None:
        """Set the list of segment IDs for the entity object"""
        if value is None:
            self.internal_segments = ""
        else:
            self.internal_segments = ",".join(str(v) for v in value)

    internal_atom_type_unique: StringProperty(  # type: ignore
        name="Unique Atom Types",
        description="A list of the unique atom types for the entity object",
        default="",
        options={"HIDDEN"},
    )

    @property
    def atom_type_unique(self) -> list[str]:
        """Return a list of the unique atom types for the entity object"""
        return (
            self.internal_atom_type_unique.split(",")
            if self.internal_atom_type_unique
            else []
        )

    @atom_type_unique.setter
    def atom_type_unique(self, value: list[str] | None) -> None:
        """Set the list of the unique atom types for the entity object"""
        if value is None:
            self.internal_atom_type_unique = ""
        else:
            self.internal_atom_type_unique = ",".join(str(v) for v in value)

    internal_categories: StringProperty(  # type: ignore
        name="Categories",
        description="Per-column category labels for the entity object",
        default="",
        options={"HIDDEN"},
    )

    @property
    def categories(self) -> dict[str, list[str]]:
        """Return the per-column category labels for the entity object"""
        if not self.internal_categories:
            return {}
        return json.loads(self.internal_categories)

    @categories.setter
    def categories(self, value: dict[str, list[str]] | str | None) -> None:
        """Set the per-column category labels for the entity object"""
        if value is None or value == "":
            self.internal_categories = ""
        elif isinstance(value, str):
            self.internal_categories = value
        else:
            self.internal_categories = json.dumps(value)

    @property
    def is_entity(self) -> bool:
        """Whether this object is a Molecular Nodes entity."""
        return self.entity_type != "None"

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

    def node(self) -> NamedAttribute:
        """Creates and returns a NamedAttribute node inside the active node tree for this selection."""
        return NamedAttribute.boolean(name=self.name)

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
    DSSPProperties,
    MolecularNodesObjectProperties,
    MolecularNodesSceneProperties,
    TrajectorySelectionItem,  # item has to be registered the ListUI and to work properly
    MN_UL_TrajectorySelectionListUI,
    MN_OT_Universe_Selection_Add,
    MN_OT_Universe_Selection_Delete,
]
