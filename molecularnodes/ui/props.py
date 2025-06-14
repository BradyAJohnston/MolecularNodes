import bpy
from bpy.props import BoolProperty, EnumProperty, IntProperty, StringProperty
from bpy.types import PropertyGroup
from ..handlers import _update_entities
from ..session import get_session
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


class MolecularNodesSceneProperties(PropertyGroup):
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

    import_format_alphafold: EnumProperty(  # type: ignore
        name="Format",
        description="Format to download as from the PDB",
        items=(
            # ("bcif", ".bcif", "Binary compressed .cif file, fastest for downloading"),
            ("cif", ".cif", "The new standard of .cif / .mmcif"),
            ("pdb", ".pdb", "The classic (and depcrecated) PDB format"),
        ),
    )

    import_format_wwpdb: EnumProperty(  # type: ignore
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
                "density_wire",
                "Wire",
                "A wire mesh surface based on the specified threshold",
                1,
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


class MolecularNodesObjectProperties(PropertyGroup):
    biological_assemblies: StringProperty(  # type: ignore
        name="Biological Assemblies",
        description="A list of biological assemblies to be created",
        default="",
    )

    entity_type: EnumProperty(  # type: ignore
        name="Entity Type",
        description="How the file was imported, dictating how MN interacts with it",
        items=(
            ("molecule", "Molecule", "A single molecule"),
            ("ensemble", "Ensemble", "A collection of molecules"),
            ("density", "Density", "An electron density map"),
            ("md", "Trajectory", "A molecular dynamics trajectory"),
            ("md-oxdna", "oxDNA Trajectory", "A oxDNA molecular dynamics trajectory "),
        ),
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

    uuid: StringProperty(  # type: ignore
        name="UUID",
        description="Unique ID for matching selection in UI to selection on python object",
        default="",
    )

    name: StringProperty(  # type: ignore
        name="Name",
        description="Name of the attribute on the mesh",
        default="custom_selection",
        update=_update_entities,
    )

    selection_str: StringProperty(  # type: ignore
        name="Selection",
        description="Selection to be applied, written in the MDAnalysis selection language",
        default="name CA",
        update=_update_entities,
    )

    updating: BoolProperty(  # type: ignore
        name="Updating",
        description="Recalculate the selection on scene frame change",
        default=True,
        # update=_selection_update_trajectories,
    )

    periodic: BoolProperty(  # type: ignore
        name="Periodic",
        description="For geometric selections, whether to account for atoms in different periodic images when searching",
        default=True,
        # update=_selection_update_trajectories,
    )

    message: StringProperty(  # type: ignore
        name="Message",
        description="Message to report back from `universe.select_atoms()`",
        default="",
    )

    immutable: BoolProperty(  # type: ignore
        name="Immutable",
        description="Whether the selection is immutable",
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
            if item.immutable:
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
        obj = context.active_object
        traj = get_session(context).match(obj)
        i = int(len(obj.mn_trajectory_selections) - 1)
        name = "selection_0"
        while True:
            if len(obj.mn_trajectory_selections) == 0:
                break
            if name in obj.mn_trajectory_selections:
                i += 1
                name = f"selection_{i}"
            else:
                break
        traj.add_selection(name=name, selection_str="all")
        obj.mn["list_index"] = i

        return {"FINISHED"}


class MN_OT_Universe_Selection_Delete(bpy.types.Operator):
    bl_idname = "mda.delete_item"
    bl_label = "-"
    bl_description = "Delete the given boolean selection from the universe"

    @classmethod
    def poll(cls, context):
        return context.active_object.mn_trajectory_selections

    def execute(self, context):
        obj = context.active_object
        index = obj.mn.trajectory_selection_index
        traj = get_session(context).match(obj)
        names = [s.name for s in obj.mn_trajectory_selections]
        traj.remove_selection(names[index])
        obj.mn.trajectory_selection_index = int(
            max(min(index, len(obj.mn_trajectory_selections) - 1), 0)
        )

        return {"FINISHED"}


CLASSES = [
    MolecularNodesObjectProperties,
    MolecularNodesSceneProperties,
    TrajectorySelectionItem,  # item has to be registered the ListUI and to work properly
    MN_UL_TrajectorySelectionListUI,
    MN_OT_Universe_Selection_Add,
    MN_OT_Universe_Selection_Delete,
]
