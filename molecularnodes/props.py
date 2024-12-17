import bpy
from bpy.types import PropertyGroup
from bpy.props import IntProperty, BoolProperty, EnumProperty, StringProperty
from .handlers import _update_trajectories
from .style import STYLE_ITEMS


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

bpy.types.Object.uuid = StringProperty(  # type: ignore
    name="UUID",
    description="Unique ID for referencing the required objects in the MNSession",
    default="",
)


class MolecularNodesSceneProperties(PropertyGroup):
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


class MolecularNodesObjectProperties(PropertyGroup):
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
        update=_update_trajectories,
        min=0,
    )
    update_with_scene: BoolProperty(  # type: ignore
        name="Update with Scene",
        description="Update the trajectory with the scene frame",
        default=True,
        update=_update_trajectories,
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
        description="Average the position this number of frames either side of the current frame",
        default=0,
        update=_update_trajectories,
        min=0,
        soft_max=5,
    )
    correct_periodic: BoolProperty(  # type: ignore
        name="Correct",
        description="Correct for periodic boundary crossing when using interpolation or averaging. Assumes cubic dimensions and only works if the unit cell is orthorhombic",
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


CLASSES = [MolecularNodesObjectProperties, MolecularNodesSceneProperties]
