import json
import math
from typing import Annotated, cast
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
from nodebpy.nodes.geometry import NamedAttribute
from ..blender.registered import object_props, object_uuid
from ..blender.utils import annotated_props, set_object_visibility
from ..entities import Density, Molecule
from ..entities.base import EntityType
from ..handlers import _update_entities, _update_selections
from ..session import get_entity, get_session

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
    ("md-oxdna", "oxDNA Trajectory", "A oxDNA molecular dynamics trajectory"),
    (
        "md-streaming",
        "Streaming Trajectory",
        "A streaming IMD molecular dynamics trajectory",
    ),
    ("ensemble-star", "Star Ensemble", "A starfile ensemble"),
    ("ensemble-cellpack", "CellPack Ensemble", "A CellPack model ensemble"),
)


def _get_frame(self):
    return self.get("frame", 0)


def _set_frame(self, frame):
    # the manual frame is stepped through subframes like the scene frame, so the
    # last frame that maps onto the trajectory is the last universe frame expanded
    # by the subframes. `_mapped_uframes` clamps the universe frame regardless;
    # this keeps the value shown in the UI on the trajectory
    last_frame = max((self.n_frames - 1) * (self.subframes + 1), 0)
    if frame > last_frame:
        frame = last_frame
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
    entity = get_session().get(object_uuid(obj))
    # ensembles have no annotations
    if isinstance(entity, (Molecule, Density)):
        entity.annotations._update_annotation_object()


def _get_entities_active_index(self) -> int:
    """
    Derive the list selection from the scene's active object.

    The active object is the single source of truth for the Entities list, so
    selecting an entity object in the outliner or viewport is mirrored in the
    list. Returns -1 when the active object is not a molecular entity.
    """
    view_layer = bpy.context.view_layer
    if view_layer is None:
        return -1
    obj = view_layer.objects.active
    if obj is None or not object_props(obj).is_entity:
        return -1
    return bpy.data.objects.find(obj.name)


def _set_entities_active_index(self, index: int) -> None:
    """set callback: mirror the list selection into the scene selection"""
    # the index points into bpy.data.objects (the Entities list data source)
    if index < 0 or index >= len(bpy.data.objects):
        return
    entity_object = bpy.data.objects[index]
    if not object_props(entity_object).is_entity:
        return
    context = bpy.context
    view_layer = context.view_layer
    if view_layer is None or entity_object.name not in view_layer.objects:
        return
    for obj in list(context.selected_objects):
        obj.select_set(False)
    view_layer.objects.active = entity_object
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
        entity.calculations._dssp.show_none()
        _update_entities(self, context)
    elif display_option == "per-frame":
        entity.calculations._dssp.show_per_frame()
        _update_entities(self, context)
    elif display_option == "sliding-window-average":
        sw_threshold = self.sw_threshold if self.apply_sw_threshold else None
        entity.calculations._dssp.show_sliding_window_average(
            self.window_size, sw_threshold
        )
        _update_entities(self, context)
    else:
        self.applied = False


def _update_dssp_applied(self, context):
    if self.applied:
        _update_entities(self, context)


@annotated_props
class DSSPProperties(bpy.types.PropertyGroup):
    selection: Annotated[
        str,
        StringProperty(
            name="Selection",
            description=(
                "MDAnalysis selection string limiting the atoms DSSP is run on. "
                "Use this to exclude termini or non-standard residues that lack "
                "the required (N, CA, C, O) atoms"
            ),
            default="protein",
        ),
    ]
    display_option: Annotated[
        str,
        EnumProperty(
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
        ),
    ]
    display_option_streaming: Annotated[
        str,
        EnumProperty(
            name="Display",
            description="Options to display secondary structures",
            items=(
                ("none", "None", "Do not show secondary structures"),
                ("per-frame", "Per Frame", "Secondary structures calculated per frame"),
            ),
            default="per-frame",
            update=_update_dssp_display_option,
        ),
    ]
    window_size: Annotated[
        int,
        IntProperty(
            name="Window Size",
            description="Number of frames in the sliding window",
            min=1,
            soft_max=10,
            default=5,
            update=_update_dssp_display_option,
        ),
    ]
    apply_sw_threshold: Annotated[
        bool,
        BoolProperty(
            name="Apply Threshold",
            description="Apply a threshold comparison to calculated mean",
            default=False,
            update=_update_dssp_display_option,
        ),
    ]
    sw_threshold: Annotated[
        float,
        FloatProperty(
            name="Threshold",
            description="Threshold fraction of frames for sliding window average",
            subtype="FACTOR",
            min=0.0,
            max=1.0,
            default=0.5,
            update=_update_dssp_display_option,
        ),
    ]
    apply_ta_threshold: Annotated[
        bool,
        BoolProperty(
            name="Apply Threshold",
            description="Apply a threshold comparison to calculated mean",
            default=False,
            update=_update_dssp_display_option,
        ),
    ]
    ta_threshold: Annotated[
        float,
        FloatProperty(
            name="Threshold",
            description="Threshold fraction of frames for trajectory average",
            subtype="FACTOR",
            min=0.0,
            max=1.0,
            default=0.5,
            update=_update_dssp_display_option,
        ),
    ]
    applied: Annotated[
        bool,
        BoolProperty(
            default=True,
            update=_update_dssp_applied,
        ),
    ]
    cancelling: Annotated[bool, BoolProperty(default=False)]


def _update_hbonds_cutoffs(self, context):
    entity = context.scene.MNSession.get(self.id_data.uuid)
    if entity is None:
        return
    entity.calculations._hbonds._update_from_props()


@annotated_props
class HBondProperties(bpy.types.PropertyGroup):
    d_a_cutoff: Annotated[
        float,
        FloatProperty(
            name="Distance",
            description="Maximum distance between the donor and acceptor in Angstroms",
            min=0.0,
            soft_max=5.0,
            default=3.0,
            update=_update_hbonds_cutoffs,
        ),
    ]
    d_h_a_angle_cutoff: Annotated[
        float,
        FloatProperty(
            name="Angle",
            description="Minimum donor-hydrogen-acceptor angle",
            subtype="ANGLE",
            min=0.0,
            max=math.pi,
            default=math.radians(150.0),
            update=_update_hbonds_cutoffs,
        ),
    ]
    d_h_cutoff: Annotated[
        float,
        FloatProperty(
            name="Donor Distance",
            description="Maximum distance between a hydrogen and its donor in Angstroms, which pairs them",
            min=0.0,
            soft_max=2.0,
            default=1.2,
            update=_update_hbonds_cutoffs,
        ),
    ]
    selection: Annotated[
        str,
        StringProperty(
            name="Selection",
            description=(
                "MDAnalysis selection string limiting the atoms considered as "
                "hydrogens, donors and acceptors"
            ),
            default="all",
        ),
    ]
    kwargs: Annotated[
        str,
        StringProperty(
            name="Arguments",
            description="Other arguments to HydrogenBondAnalysis than the cutoffs, as JSON, restored with the file",
            default="{}",
            options={"HIDDEN"},
        ),
    ]


FRAME_ATTRIBUTE_SOURCES = (
    (
        "TRAJECTORY",
        "Trajectory",
        "Read from each frame of the trajectory",
        "FILE_MOVIE",
        0,
    ),
    ("SELECTION", "Selection", "A selection of the molecule", "RESTRICT_SELECT_OFF", 1),
    ("CALCULATION", "Calculation", "Calculated from the universe", "DRIVER", 2),
)


@annotated_props
class FrameAttributeItem(bpy.types.PropertyGroup):
    """An attribute that is written on frame changes, and where it comes from."""

    __slots__ = []

    name: Annotated[
        str,
        StringProperty(
            name="Attribute Name",
            description="Name of the attribute on the mesh",
        ),
    ]
    source: Annotated[
        str,
        EnumProperty(
            name="Source",
            description="Where the values of the attribute come from",
            items=FRAME_ATTRIBUTE_SOURCES,
            default="TRAJECTORY",
        ),
    ]


@annotated_props
class CalculationItem(bpy.types.PropertyGroup):
    """A calculation that writes an attribute from the universe."""

    __slots__ = []

    name: Annotated[
        str,
        StringProperty(
            name="Attribute Name",
            description="Name of the attribute the calculation writes",
        ),
    ]
    updating: Annotated[
        bool,
        BoolProperty(
            name="Updating",
            description="Recalculate on every frame change, rather than keeping the values from when it last ran",
            default=True,
            update=_update_entities,
        ),
    ]
    object: Annotated[
        bpy.types.Object | None,
        PointerProperty(
            type=bpy.types.Object,
            name="Object",
            description="Object the calculation writes to, for one that writes an object rather than an attribute",
        ),
    ]


@annotated_props
class MolecularNodesSceneProperties(bpy.types.PropertyGroup):
    __slots__ = []
    entities_active_index: Annotated[
        int,
        IntProperty(
            name="Active entity index",
            description="Index into bpy.data.objects of the entity object active in"
            " the Entities list, derived from the scene's active object",
            default=-1,
            get=_get_entities_active_index,
            set=_set_entities_active_index,
        ),
    ]

    is_updating: Annotated[
        bool,
        BoolProperty(
            name="Updating",
            description="Currently updating data in the scene, don't trigger more updates",
            default=False,
        ),
    ]


def _update_annotations_visibility(self, context):
    entity = context.scene.MNSession.get(self.id_data.uuid)
    if entity is not None:
        if self.annotations_visible:
            entity.annotations._draw_handler_add()
        else:
            entity.annotations._draw_handler_remove()
        entity.annotations._update_annotation_object()


@annotated_props
class MolecularNodesObjectProperties(bpy.types.PropertyGroup):
    __slots__ = []
    styles_active_index: Annotated[int, IntProperty(default=-1)]
    annotations_active_index: Annotated[int, IntProperty(default=-1)]
    annotations_next_index: Annotated[int, IntProperty(default=0)]

    annotations_visible: Annotated[
        bool,
        BoolProperty(
            name="Visible",
            description="Visibility of all annotations",
            default=True,
            update=_update_annotations_visibility,
        ),
    ]

    internal_biological_assemblies: Annotated[
        str,
        StringProperty(
            name="Biological Assemblies",
            description="A list of biological assemblies to be created",
            default="",
            options={"HIDDEN"},
        ),
    ]

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

    entity_type: Annotated[
        str,
        EnumProperty(
            name="Entity Type",
            description="How the file was imported, dictating how MN interacts with it",
            items=ENTITY_ITEMS,
            default="None",
        ),
    ]

    visible: Annotated[
        bool,
        BoolProperty(
            name="Visible",
            description="Visibility of the entity object in the viewport and renders",
            get=_get_object_visibility,
            set=_set_object_visibility,
        ),
    ]

    code: Annotated[
        str,
        StringProperty(
            name="PDB",
            description="PDB code used to download this structure",
            maxlen=12,
            options={"HIDDEN"},
        ),
    ]
    filepath: Annotated[
        str,
        StringProperty(
            name="File Path",
            description="Source file this entity was loaded from, used to reload it",
            subtype="FILE_PATH",
            default="",
            options={"HIDDEN"},
        ),
    ]
    data_block: Annotated[
        str,
        StringProperty(
            name="Data Block",
            description="Data block of a multi-structure CIF file this entity was loaded from",
            default="",
            options={"HIDDEN"},
        ),
    ]
    database: Annotated[
        str,
        StringProperty(
            name="Database",
            description="Database this structure was fetched from, used to reload it",
            default="",
            options={"HIDDEN"},
        ),
    ]

    internal_chain_ids: Annotated[
        str,
        StringProperty(
            name="Chain IDs",
            description="A list of chain IDs for the entity object",
            default="",
            options={"HIDDEN"},
        ),
    ]

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

    internal_entity_ids: Annotated[
        str,
        StringProperty(
            name="Entity IDs",
            description="A list of entity IDs for the entity object",
            default="",
            options={"HIDDEN"},
        ),
    ]

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

    internal_segments: Annotated[
        str,
        StringProperty(
            name="Segments",
            description="A list of segment IDs for the entity object",
            default="",
            options={"HIDDEN"},
        ),
    ]

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

    internal_atom_type_unique: Annotated[
        str,
        StringProperty(
            name="Unique Atom Types",
            description="A list of the unique atom types for the entity object",
            default="",
            options={"HIDDEN"},
        ),
    ]

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

    internal_categories: Annotated[
        str,
        StringProperty(
            name="Categories",
            description="Per-column category labels for the entity object",
            default="",
            options={"HIDDEN"},
        ),
    ]

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

    trajectory_selection_index: Annotated[
        int,
        IntProperty(
            name="Index of selection",
            description="Index of selection, that is selected for the UI",
            default=0,
        ),
    ]
    frame_hidden: Annotated[
        int,
        IntProperty(
            name="Frame",
            description="Frame of the loaded trajectory",
            default=0,
            min=0,
        ),
    ]
    frame: Annotated[
        int,
        IntProperty(
            name="Frame",
            description="Manual frame of the loaded trajectory. Stepped through subframes like the scene frame, with the resolved trajectory frame shown alongside",
            default=0,
            min=0,
            set=_set_frame,
            get=_get_frame,
        ),
    ]
    n_frames: Annotated[
        int,
        IntProperty(
            name="Number of Frames",
            description="Number of frames in the loaded trajectory",
            default=0,
            min=0,
        ),
    ]
    update_with_scene: Annotated[
        bool,
        BoolProperty(
            name="Update with Scene",
            description="Update the trajectory with the scene frame",
            default=True,
            update=_update_entities,
        ),
    ]
    subframes: Annotated[
        int,
        IntProperty(
            name="Subframes",
            description="Number of subframes to insert between frames of the loaded trajectory",
            default=0,
            update=_update_entities,
            min=0,
        ),
    ]
    offset: Annotated[
        int,
        IntProperty(
            name="Offset",
            description="Offset the starting playback for the trajectory on the timeline. Positive starts the playback later than frame 0, negative starts it earlier than frame 0",
            default=0,
            update=_update_entities,
        ),
    ]

    interpolate: Annotated[
        bool,
        BoolProperty(
            name="Interpolate",
            description="Whether to interpolate when using subframes",
            default=True,
            update=_update_entities,
        ),
    ]
    average: Annotated[
        int,
        IntProperty(
            name="Average",
            description="Average the positions across this number of frames on either side of the current frame",
            default=0,
            update=_update_entities,
            min=0,
            soft_max=4,
        ),
    ]
    correct_periodic: Annotated[
        bool,
        BoolProperty(
            name="Correct",
            description="Correct for periodic boundary crossing when using interpolation or averaging, moving atoms to the nearest periodic image. Works with orthorhombic and triclinic unit cells",
            default=False,
            update=_update_entities,
        ),
    ]
    frame_attributes: Annotated[
        "bpy.types.bpy_prop_collection_idprop[FrameAttributeItem]",
        CollectionProperty(
            type=FrameAttributeItem,
            name="Per-frame Attributes",
            description="The attributes written on frame changes, from the trajectory, selections and calculations",
            options={"HIDDEN"},
        ),
    ]
    calculations: Annotated[
        "bpy.types.bpy_prop_collection_idprop[CalculationItem]",
        CollectionProperty(
            type=CalculationItem,
            name="Calculations",
            description="Calculations that write attributes from the universe, and whether they update on frame changes",
            options={"HIDDEN"},
        ),
    ]
    frame_attributes_index: Annotated[
        int,
        IntProperty(
            name="Active Per-frame Attribute",
            default=0,
            options={"HIDDEN"},
        ),
    ]
    frame_attributes_filter: Annotated[
        str,
        EnumProperty(
            name="Show",
            description="Which per-frame attributes to list, by where they come from",
            items=(("ALL", "All", "List every per-frame attribute", "NONE", 3),)
            + FRAME_ATTRIBUTE_SOURCES,
            default="ALL",
            options={"HIDDEN"},
        ),
    ]
    filepath_trajectory: Annotated[
        str,
        StringProperty(
            name="Trajectory",
            description="Filepath for the `trajectory` of the Object",
            subtype="FILE_PATH",
            default="",
        ),
    ]
    filepath_topology: Annotated[
        str,
        StringProperty(
            name="Topology",
            description="Filepath for the Topology of the Object",
            subtype="FILE_PATH",
            default="",
        ),
    ]
    dssp: Annotated[DSSPProperties, PointerProperty(type=DSSPProperties)]
    hbonds: Annotated[HBondProperties, PointerProperty(type=HBondProperties)]


@annotated_props
class TrajectorySelectionItem(bpy.types.PropertyGroup):
    """Group of properties for custom selections for MDAnalysis import."""

    __slots__ = []

    def node(self) -> NamedAttribute:
        """Creates and returns a NamedAttribute node inside the active node tree for this selection."""
        return NamedAttribute.boolean(name=self.name)

    name: Annotated[
        str,
        StringProperty(
            name="Attribute Name",
            description="Name of the attribute that will be created when storing on the mesh",
        ),
    ]

    string: Annotated[
        str,
        StringProperty(
            name="Selection",
            description="Selection to be applied, written in the MDAnalysis selection language",
            default="name CA",
            update=_update_selections,
        ),
    ]

    previous_string: Annotated[str, StringProperty()]
    previous_updating: Annotated[bool, BoolProperty(default=True)]
    previous_periodic: Annotated[bool, BoolProperty(default=True)]

    updating: Annotated[
        bool,
        BoolProperty(
            name="Updating",
            description="Potential recalculate the selection when the scene frame changes",
            default=True,
            update=_update_selections,
        ),
    ]

    periodic: Annotated[
        bool,
        BoolProperty(
            name="Periodic",
            description="For geometric selections, whether to account for atoms in different periodic images when searching",
            default=True,
            update=_update_selections,
        ),
    ]

    message: Annotated[
        str,
        StringProperty(
            name="Message",
            description="Message to report back from `universe.select_atoms()`",
            default="",
        ),
    ]

    from_atomgroup: Annotated[
        bool,
        BoolProperty(
            name="From AtomGroup",
            description="If the UI item has been created from an existing AtomGroup. Will prevent editing in the UI by a user.",
            default=False,
        ),
    ]


@annotated_props
class MN_OT_Universe_Selection_Add(bpy.types.Operator):
    "Add a new custom selection to a trajectory"

    bl_idname = "mn.trajectory_selection_add"
    bl_label = "+"
    bl_description = "Add a new boolean attribute for the given MDA selection string"

    uuid: Annotated[str, StringProperty()]

    def execute(self, context):
        if self.uuid:
            traj = context.scene.MNSession.get(self.uuid)
        else:
            traj = get_entity(context)
        # selections are only added from the panels of trajectories
        traj = cast(Molecule, traj)
        item = traj.selections.from_string("all")
        traj.selections.ui_index = max(0, len(traj.selections) - 1)
        # make the new selection the active row of the Attributes list, so its
        # selection string shows for editing
        mn = traj.props
        if mn.frame_attributes_filter not in ("ALL", "SELECTION"):
            mn.frame_attributes_filter = "ALL"
        index = mn.frame_attributes.find(item.name)
        if index >= 0:
            mn.frame_attributes_index = index
        return {"FINISHED"}


class MN_OT_Universe_Selection_Delete(bpy.types.Operator):
    bl_idname = "mn.trajectory_selection_remove"
    bl_label = "-"
    bl_description = "Delete the given boolean selection from the universe"

    @classmethod
    def poll(cls, context):
        return context.active_object.mn_trajectory_selections

    def execute(self, context):
        traj = cast(Molecule, get_entity(context))
        index = traj.selections.ui_index
        traj.selections.remove(index)

        # the length of items in the list has changed, set the currently selected index
        # to a new value. Ensure it is between 0 and the length of the items in the list
        traj.selections.ui_index = max(0, min(index, len(traj.selections.ui_items) - 1))

        return {"FINISHED"}


CLASSES = [
    DSSPProperties,
    HBondProperties,
    FrameAttributeItem,
    CalculationItem,
    MolecularNodesObjectProperties,
    MolecularNodesSceneProperties,
    TrajectorySelectionItem,  # item has to be registered the ListUI and to work properly
    MN_OT_Universe_Selection_Add,
    MN_OT_Universe_Selection_Delete,
]
