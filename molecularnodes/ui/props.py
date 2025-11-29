import bpy
from bpy.props import (  # type: ignore
    BoolProperty,
    CollectionProperty,
    EnumProperty,
    IntProperty,
    StringProperty,
)
from bpy.types import PropertyGroup  # type: ignore
from databpy.object import LinkedObjectError
from ..blender.utils import set_object_visibility
from ..handlers import _update_entities
from ..session import get_entity
from .style import STYLE_ITEMS

uuid_property = StringProperty(  # type: ignore
    name="UUID",
    description="Unique ID for referencing the required objects in the MNSession",
    default="",
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
    visible: BoolProperty(
        name="visible",
        description="Visibility of the entity",
        default=True,
        get=_get_entity_visibility,
        set=_set_entity_visibility,
    )  # type: ignore


class MolecularNodesSceneProperties(PropertyGroup):
    __slots__ = []
    entities: CollectionProperty(name="Entities", type=EntityProperties)  # type: ignore
    entities_active_index: IntProperty(
        name="Active entity index",
        default=-1,
        update=_entities_active_index_callback,
    )  # type: ignore

    import_del_hydrogen: BoolProperty(  # type: ignore
        name="Remove Hydrogens",
        description="Remove the hydrogens from a structure on import",
        default=False,
    )

    import_local_path: StringProperty(  # type: ignore
        name="File",
        description="File path of the structure to open",
        options={"TEXTEDIT_UPDATE"},
        subtype="FILE_PATH",
        maxlen=0,
    )

    import_code_alphafold: StringProperty(  # type: ignore
        name="UniProt ID",
        description="The UniProt ID to use for downloading from the AlphaFold databse",
        options={"TEXTEDIT_UPDATE"},
    )

    import_format_fetch: EnumProperty(  # type: ignore
        name="Format",
        description="Format to download as from the PDB",
        items=(
            ("bcif", ".bcif", "Binary compressed .cif file, fastest for downloading"),
            ("cif", ".cif", "The new standard of .cif / .mmcif"),
            ("pdb", ".pdb", "The classic (and depcrecated) PDB format"),
        ),
    )

    import_code_pdb: StringProperty(  # type: ignore
        name="PDB",
        description="The PDB code to download and import",
        options={"TEXTEDIT_UPDATE"},
        maxlen=4,
    )

    is_updating: BoolProperty(  # type: ignore
        name="Updating",
        description="Currently updating data in the scene, don't trigger more updates",
        default=False,
    )

    import_centre: BoolProperty(  # type: ignore
        name="Centre Structure",
        description="Move the imported Molecule on the World Origin",
        default=False,
    )

    import_centre_type: EnumProperty(  # type: ignore
        name="Method",
        default="mass",
        items=(
            (
                "mass",
                "Mass",
                "Adjust the structure's centre of mass to be at the world origin",
                1,
            ),
            (
                "centroid",
                "Centroid",
                "Adjust the structure's centroid (centre of geometry) to be at the world origin",
                2,
            ),
        ),
    )

    import_node_setup: BoolProperty(  # type: ignore
        name="Setup Nodes",
        default=True,
        description="Create and set up a Geometry Nodes tree on import",
    )

    import_build_assembly: BoolProperty(  # type: ignore
        name="Build Assembly",
        description="Add a node to build the biological assembly on import",
        default=False,
    )

    import_remove_solvent: BoolProperty(  # type: ignore
        name="Remove Solvent",
        description="Delete the solvent from the structure on import",
        default=True,
    )

    import_oxdna_topology: StringProperty(  # type: ignore
        name="Toplogy",
        description="File path for the topology to import (.top)",
        subtype="FILE_PATH",
    )
    import_oxdna_trajectory: StringProperty(  # type: ignore
        name="Trajectory",
        description="File path for the trajectory to import (.oxdna / .dat)",
        subtype="FILE_PATH",
    )
    import_oxdna_name: StringProperty(  # type: ignore
        name="Name", description="Name of the created object.", default="NewOrigami"
    )
    import_style: EnumProperty(  # type: ignore
        name="Style",
        description="Default style for importing",
        items=STYLE_ITEMS,
        default="spheres",
    )
    import_md_topology: StringProperty(  # type: ignore
        name="Topology",
        description="File path for the toplogy file for the trajectory",
        subtype="FILE_PATH",
        maxlen=0,
    )
    import_md_trajectory: StringProperty(  # type: ignore
        name="Trajectory",
        description="File path for the trajectory file for the trajectory",
        subtype="FILE_PATH",
        maxlen=0,
    )
    import_md_name: StringProperty(  # type: ignore
        name="Name",
        description="Name of the molecule on import",
        default="NewTrajectory",
        maxlen=0,
    )
    import_density_invert: BoolProperty(  # type: ignore
        name="Invert Data",
        description="Invert the values in the map. Low becomes high, high becomes low.",
        default=False,
    )
    import_density_center: BoolProperty(  # type: ignore
        name="Center Density",
        description="Translate the density so that the center of the box is at the origin.",
        default=False,
    )
    import_density_overwrite: BoolProperty(  # type: ignore
        name="Overwrite Intermediate File",
        description="Overwrite generated intermediate .vdb file.",
        default=False,
    )
    import_density: StringProperty(  # type: ignore
        name="File",
        description="File path for the map file.",
        subtype="FILE_PATH",
        maxlen=0,
    )

    import_density_style: EnumProperty(  # type: ignore
        name="Style",
        items=(
            (
                "density_surface",
                "Surface",
                "A mesh surface based on the specified threshold",
                0,
            ),
            (
                "density_iso_surface",
                "ISO Surface",
                "A mesh surface based on the specified iso value",
                1,
            ),
            (
                "density_wire",
                "Wire",
                "A wire mesh surface based on the specified threshold",
                2,
            ),
        ),
    )

    panel_selection: bpy.props.EnumProperty(  # type: ignore
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

    panel_import_type: bpy.props.EnumProperty(  # type: ignore
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
    import_star_file_path: StringProperty(  # type: ignore
        name="File",
        description="File path for the `.star` file to import.",
        subtype="FILE_PATH",
        maxlen=0,
    )
    import_cell_pack_path: bpy.props.StringProperty(  # type: ignore
        name="File",
        description="File to import (.cif, .bcif)",
        subtype="FILE_PATH",
        maxlen=0,
    )


def _update_annotations_visibility(self, context):
    entity = context.scene.MNSession.get(self.id_data.uuid)
    if entity is not None:
        entity.annotations._update_annotation_object()


class MolecularNodesObjectProperties(PropertyGroup):
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
        items=(
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
        ),
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
    MolecularNodesObjectProperties,
    MolecularNodesSceneProperties,
    TrajectorySelectionItem,  # item has to be registered the ListUI and to work properly
    MN_UL_TrajectorySelectionListUI,
    MN_OT_Universe_Selection_Add,
    MN_OT_Universe_Selection_Delete,
]
