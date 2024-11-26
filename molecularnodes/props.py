import bpy
from bpy.props import IntProperty, BoolProperty, EnumProperty, StringProperty
from .entities.trajectory.handlers import _update_trajectories

bpy.types.Scene.MN_import_centre = BoolProperty(
    name="Centre Structure",
    description="Move the imported Molecule on the World Origin",
    default=False,
)

bpy.types.Scene.MN_centre_type = EnumProperty(
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

bpy.types.Scene.MN_import_del_solvent = BoolProperty(
    name="Remove Solvent",
    description="Delete the solvent from the structure on import",
    default=True,
)
bpy.types.Scene.MN_import_panel_selection = IntProperty(
    name="MN_import_panel_selection",
    description="Import Panel Selection",
    subtype="NONE",
    default=0,
)
bpy.types.Scene.MN_import_build_assembly = BoolProperty(
    name="Build Assembly", default=False
)
bpy.types.Scene.MN_import_node_setup = BoolProperty(
    name="Setup Nodes",
    default=True,
    description="Create and set up a Geometry Nodes tree on import",
)


class MolecularNodesObjectProperties(bpy.types.PropertyGroup):
    molecule_type: StringProperty(  # type: ignore
        name="Molecular Type",
        description="How the file was imported, dictating how MN interacts with it",
        default="",
    )
    uuid: StringProperty(  # type: ignore
        name="UUID",
        description="Unique ID for referencing the required objects in the MNSession",
        default="",
    )
    pdb_code: StringProperty(  # type: ignore
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
    subframes: IntProperty(  # type: ignore
        name="Subframes",
        description="Number of subframes to insert between frames of the loaded trajectory",
        default=0,
        update=_update_trajectories,
        min=0,
    )
    offset: IntProperty(  # type: ignore
        name="Offset",
        description="Offset the starting playback for the trajectory on the timeine. Positive starts the playback later than frame 0, negative starts it earlier than frame 0",
        default=0,
        update=_update_trajectories,
    )
    interpolate: BoolProperty(  # type: ignore
        name="Interpolate",
        description="Whether to interpolate when using subframes",
        default=True,
        update=_update_trajectories,
    )
    average: IntProperty(  # type: ignore
        name="Average",
        description="Average values between frame +/- the number of selected frames",
        default=0,
        update=_update_trajectories,
        min=0,
        soft_max=5,
    )
    correct_periodic: BoolProperty(  # type: ignore
        name="Correct",
        description="Correct for periodic boundary crossing when using interpolation. Assumes cubic dimensions",
        default=False,
        update=_update_trajectories,
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
