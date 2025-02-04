import bpy
from bpy.types import PropertyGroup
from bpy.props import IntProperty, BoolProperty, EnumProperty, StringProperty
from .handlers import _selection_update_trajectories, _udpate_entities
from .style import STYLE_ITEMS


uuid_property = StringProperty(  # type: ignore
    name="UUID",
    description="Unique ID for referencing the required objects in the MNSession",
    default="",
)


class MolecularNodesSceneProperties(PropertyGroup):
    import_centre: BoolProperty(  # type: ignore
        name="Centre Structure",
        description="Move the imported Molecule on the World Origin",
        default=False,
    )

    centre_type: EnumProperty(  # type: ignore
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

    import_del_solvent: BoolProperty(  # type: ignore
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
        update=_udpate_entities,
        min=0,
    )
    update_with_scene: BoolProperty(  # type: ignore
        name="Update with Scene",
        description="Update the trajectory with the scene frame",
        default=True,
        update=_udpate_entities,
    )
    subframes: IntProperty(  # type: ignore
        name="Subframes",
        description="Number of subframes to insert between frames of the loaded trajectory",
        default=0,
        update=_udpate_entities,
        min=0,
    )
    offset: IntProperty(  # type: ignore
        name="Offset",
        description="Offset the starting playback for the trajectory on the timeine. Positive starts the playback later than frame 0, negative starts it earlier than frame 0",
        default=0,
        update=_udpate_entities,
    )
    interpolate: BoolProperty(  # type: ignore
        name="Interpolate",
        description="Whether to interpolate when using subframes",
        default=True,
        update=_udpate_entities,
    )
    average: IntProperty(  # type: ignore
        name="Average",
        description="Average the position this number of frames either side of the current frame",
        default=0,
        update=_udpate_entities,
        min=0,
        soft_max=5,
    )
    correct_periodic: BoolProperty(  # type: ignore
        name="Correct",
        description="Correct for periodic boundary crossing when using interpolation or averaging. Assumes cubic dimensions and only works if the unit cell is orthorhombic",
        default=False,
        update=_udpate_entities,
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

    name: StringProperty(  # type: ignore
        name="Name",
        description="Name of the attribute on the mesh",
        default="custom_selection",
        update=_selection_update_trajectories,
    )

    selection_str: StringProperty(  # type: ignore
        name="Selection",
        description="Selection to be applied, written in the MDAnalysis selection language",
        default="name CA",
        update=_selection_update_trajectories,
    )

    updating: BoolProperty(  # type: ignore
        name="Updating",
        description="Recalculate the selection on scene frame change",
        default=True,
        update=_selection_update_trajectories,
    )

    periodic: BoolProperty(  # type: ignore
        name="Periodic",
        description="For geometric selections, whether to account for atoms in different periodic images when searching",
        default=True,
        update=_selection_update_trajectories,
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

            row.prop(item, "name", text="", emboss=False)
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
        obj.mn_trajectory_selections.add()
        i = int(len(obj.mn_trajectory_selections) - 1)
        obj.mn_trajectory_selections[i].name = f"selection_{i + 1}"
        obj.mn["list_index"] = i
        _udpate_entities(self, context)

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

        sel_list = obj.mn_trajectory_selections
        sel_list.remove(index)
        obj.mn.trajectory_selection_index = len(sel_list) - 1
        _udpate_entities(self, context)

        return {"FINISHED"}


CLASSES = [
    MolecularNodesObjectProperties,
    MolecularNodesSceneProperties,
    TrajectorySelectionItem,  # item has to be registered the ListUI and to work properly
    MN_UL_TrajectorySelectionListUI,
    MN_OT_Universe_Selection_Add,
    MN_OT_Universe_Selection_Delete,
]
