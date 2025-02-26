import bpy
from bpy.props import EnumProperty, BoolProperty, StringProperty, CollectionProperty
from pathlib import Path
import MDAnalysis as mda


from .style import STYLE_ITEMS
from ..download import FileDownloadPDBError, CACHE_DIR
from ..entities import density, ensemble, trajectory, molecule


class Import_Molecule(bpy.types.Operator):
    style: EnumProperty(  # type: ignore
        name="Style",
        default="spheres",
        description="Starting style for the structure on import",
        items=STYLE_ITEMS,
    )

    node_setup: BoolProperty(  # type: ignore
        name="Node Setup",
        default=True,
        description="Whether to setup the starting default node tree on import",
    )

    centre: BoolProperty(  # type: ignore
        name="Centre",
        description="Whether to centre the structure on import",
        default=False,
    )
    centre_type: EnumProperty(  # type: ignore
        name="Centre",
        description="Centre the structure at the world origin using the given method",
        default="mass",
        items=(
            (
                "mass",
                "Mass",
                "Adjust the structure's centre of mass to be at the world origin",
                2,
            ),
            (
                "centroid",
                "Centroid",
                "Adjust the structure's centroid (centre of geometry) to be at the world origin",
                3,
            ),
        ),
    )
    del_solvent: BoolProperty(  # type: ignore
        default=True,
        name="Delete Solvent",
        description="Remove solvent atoms from the structure on import",
    )
    assembly: BoolProperty(  # type: ignore
        default=False,
        name="Build Biological Assembly",
        description="Build the biological assembly for the structure on import",
    )

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.prop(self, "node_setup", text="")
        col = row.column()
        col.prop(self, "style")
        col.enabled = self.node_setup
        # row = layout.row()
        layout.prop(self, "centre")
        layout.prop(self, "del_solvent")
        layout.prop(self, "assembly")

        return layout


class MN_OT_Import_Molecule(Import_Molecule):
    bl_idname = "mn.import_molecule"
    bl_label = "Import a Molecule"

    directory: StringProperty(  # type: ignore
        subtype="FILE_PATH", options={"SKIP_SAVE", "HIDDEN"}
    )
    files: CollectionProperty(  # type: ignore
        type=bpy.types.OperatorFileListElement, options={"SKIP_SAVE", "HIDDEN"}
    )

    def draw(self, context):
        layout = self.layout
        layout.label(text=f"Importing {len(self.files)} molecules")
        layout = super().draw(context)

    def execute(self, context):
        if not self.directory:
            return {"CANCELLED"}

        if not self.node_setup:
            style = None
        else:
            style = self.style

        for file in self.files:
            try:
                mol = parse(os.path.join(self.directory, file.name))
                mol.create_object(
                    name=file.name,
                    centre=self.centre,
                    style=style,
                    del_solvent=self.del_solvent,
                    build_assembly=self.assembly,
                )
            except Exception as e:
                print(f"Failed importing {file}: {e}")

        return {"FINISHED"}

    def invoke(self, context, event):
        if context.area and context.area.type == "VIEW_3D":
            context.window_manager.invoke_props_dialog(self)
        else:
            context.window_manager.fileselect_add(self)
        return {"RUNNING_MODAL"}


class MN_FH_Import_Molecule(bpy.types.FileHandler):
    bl_idname = "MN_FH_import_molecule"
    bl_label = "File handler for import molecular data files."
    bl_import_operator = "mn.import_molecule"
    bl_file_extensions = ".pdb;.cif;.mmcif;.bcif;.pdbx"

    @classmethod
    def poll_drop(cls, context):
        return context.area and context.area.type == "VIEW_3D"


DOWNLOAD_FORMATS = (
    ("bcif", ".bcif", "Binary compressed .cif file, fastest for downloading"),
    ("cif", ".cif", "The new standard of .cif / .mmcif"),
    ("pdb", ".pdb", "The classic (and depcrecated) PDB format"),
)


# operator that is called by the 'button' press which calls the fetch function


class MN_OT_Import_Fetch(bpy.types.Operator):
    bl_idname = "mn.import_fetch"
    bl_label = "Fetch"
    bl_description = "Download and open a structure from the Protein Data Bank"
    bl_options = {"REGISTER", "UNDO"}

    code: StringProperty(  # type: ignore
        name="PDB",
        description="The 4-character PDB code to download",
        options={"TEXTEDIT_UPDATE"},
        maxlen=4,
    )
    file_format: EnumProperty(  # type: ignore
        name="Format",
        description="Format to download as from the PDB",
        default="bcif",
        items=DOWNLOAD_FORMATS,
    )
    node_setup: BoolProperty(  # type: ignore
        name="Setup Nodes",
        default=True,
        description="Create and set up a Geometry Nodes tree on import",
    )
    assembly: BoolProperty(  # type: ignore
        name="Build Assembly",
        description="Add a node to build the biological assembly on import",
        default=False,
    )
    style: EnumProperty(  # type: ignore
        name="Style",
        description="Default style for importing",
        items=STYLE_ITEMS,
        default="spheres",
    )
    cache_dir: StringProperty(  # type: ignore
        name="Cache Directory",
        description="Where to store the structures downloaded from the Protein Data Bank",
        default=CACHE_DIR,
        subtype="DIR_PATH",
    )
    del_solvent: BoolProperty(  # type: ignore
        name="Remove Solvent",
        description="Delete the solvent from the structure on import",
        default=True,
    )
    del_hydrogen: BoolProperty(  # type: ignore
        name="Remove Hydrogens",
        description="Remove the hydrogens from a structure on import",
        default=False,
    )

    centre: BoolProperty(  # type: ignore
        name="Centre",
        description="Centre the structure on the world origin",
        default=False,
    )

    database: EnumProperty(  # type: ignore
        name="Method",
        default="wwpdb",
        items=(
            (
                "wwpdb",
                "wwPDB",
                "The world-wide Protein Data Bank (wwPDB)",
            ),
            (
                "alphafold",
                "AlphaFold",
                "The AlphaFold computational structure database",
            ),
        ),
    )

    centre_type: EnumProperty(  # type: ignore
        name="Method",
        default="mass",
        items=(
            (
                "mass",
                "Mass",
                "Adjust the structure's centre of mass to be at the world origin",
            ),
            (
                "centroid",
                "Centroid",
                "Adjust the structure's centroid (centre of geometry) to be at the world origin",
            ),
        ),
    )

    # def invoke(self, context, event):

    #     context.window_manager.invoke_props_dialog(self)
    #     return {"RUNNING_MODAL"}

    def execute(self, context):
        try:
            mol = molecule.fetch(
                code=self.code,
                centre=self.centre_type if self.centre else None,
                del_solvent=self.del_solvent,
                del_hydrogen=self.del_hydrogen,
                style=self.style if self.node_setup else None,
                cache_dir=self.cache_dir,
                build_assembly=self.assembly,
                format=self.file_format,
            )
        except FileDownloadPDBError as e:
            self.report({"ERROR"}, str(e))
            if self.file_format == "pdb":
                self.report(
                    {"ERROR"},
                    "There may not be a `.pdb` formatted file available - try a different download format.",
                )
            return {"CANCELLED"}

        bpy.context.view_layer.objects.active = mol.object
        self.report({"INFO"}, message=f"Imported '{self.code}' as {mol.name}")

        return {"FINISHED"}


class MN_OT_Import_Protein_Local(Import_Molecule):
    bl_idname = "mn.import_local"
    bl_label = "Local"
    bl_description = "Open a local structure file"
    bl_options = {"REGISTER", "UNDO"}

    filepath: StringProperty(  # type: ignore
        name="File",
        description="File to import",
        subtype="FILE_PATH",
    )

    def execute(self, context):
        mol = molecule.parse(self.filepath)
        mol.create_object(
            name=Path(self.filepath).stem,
            style=self.style if self.node_setup else None,
            build_assembly=self.assembly,
            centre=self.centre_type if self.centre else None,
            del_solvent=self.del_solvent,
        )

        # return the good news!
        bpy.context.view_layer.objects.active = mol.object
        self.report({"INFO"}, message=f"Imported '{self.filepath}' as {mol.name}")
        return {"FINISHED"}

    def invoke(self, context, event):
        return self.execute(context)


class ImportEnsemble(bpy.types.Operator):
    filepath: StringProperty(  # type: ignore
        name="File",
        description="File path for the `.star` file to import.",
        subtype="FILE_PATH",
        maxlen=0,
    )
    node_setup: BoolProperty(  # type: ignore
        name="Setup Nodes",
        default=True,
        description="Create and set up a Geometry Nodes tree on import",
    )


class MN_OT_Import_Star_File(ImportEnsemble):
    bl_idname = "mn.import_star_file"
    bl_label = "Load"
    bl_description = (
        "Will import the given file, setting up the points to instance an object."
    )
    bl_options = {"REGISTER"}

    @classmethod
    def poll(cls, context):
        return True

    def execute(self, context):
        ensemble.load_starfile(
            file_path=self.filepath,
            node_setup=self.node_setup,
        )
        return {"FINISHED"}


class MN_OT_Import_Cell_Pack(ImportEnsemble):
    bl_idname = "mol.import_cell_pack"
    bl_label = "Load"
    bl_description = ""
    bl_options = {"REGISTER"}

    def execute(self, context):
        ens.load_cellpack(
            file_path=self.filepath,
            name=Path(self.filepath).name,
            node_setup=self.node_setup,
        )
        return {"FINISHED"}


class MN_OT_Import_Map(bpy.types.Operator):
    bl_idname = "mn.import_density"
    bl_label = "Load"
    bl_description = "Import a EM density map into Blender"
    bl_options = {"REGISTER"}

    def execute(self, context):
        scene = context.scene
        density.load(
            file_path=scene.mn.import_density,
            invert=scene.mn.import_density_invert,
            setup_nodes=scene.mn.import_node_setup,
            style=scene.mn.import_density_style,
            center=scene.mn.import_density_center,
        )
        return {"FINISHED"}


class TrajectoryImportOperator(bpy.types.Operator):
    bl_label = "Import"
    bl_description = "Will import the given file and toplogy."
    bl_options = {"REGISTER"}

    topology: StringProperty(  # type: ignore
        name="Toplogy",
        description="File path for the topology file",
        subtype="FILE_PATH",
        maxlen=0,
    )
    trajectory: StringProperty(  # type: ignore
        name="Trajectory",
        description="File path for the trajectory file",
        subtype="FILE_PATH",
        maxlen=0,
    )
    name: StringProperty(  # type: ignore
        name="Name",
        description="Name for the object that will be created and linked to the trajectory",
        default="NewOrigami",
        maxlen=0,
    )
    style: EnumProperty(  # type: ignore
        name="Style",
        description="Starting style on import",
        default="spheres",
        items=STYLE_ITEMS,
    )


class MN_OT_Reload_Trajectory(bpy.types.Operator):
    bl_idname = "mn.reload_trajectory"
    bl_label = "Reload Trajectory"
    bl_description = (
        "Reload the `mda.UNiverse` of the current Object to renable updating"
    )
    bl_options = {"REGISTER", "UNDO"}

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        traj = context.scene.MNSession.match(obj)
        return not traj

    def execute(self, context):
        obj = context.active_object
        session: MNSession = context.scene.MNSession
        topo = obj.mn.filepath_topology
        traj = obj.mn.filepath_trajectory

        if "oxdna" in obj.mn.entity_type:
            uni = mda.Universe(
                topo,
                traj,
                topology_format=trajectory.oxdna.OXDNAParser,
                format=trajectory.oxdna.OXDNAReader,
            )
            traj = trajectory.oxdna.OXDNA(uni)
        else:
            traj = trajectory.load(topo, traj)

        traj.object = obj
        traj.set_frame(context.scene.frame_current)
        return {"FINISHED"}


class MN_OT_Import_Trajectory(bpy.types.Operator):
    bl_idname = "mn.import_trajectory"
    bl_label = "Import Protein MD"
    bl_description = "Load molecular dynamics trajectory"
    bl_options = {"REGISTER", "UNDO"}

    topology: StringProperty(  # type: ignore
        name="Topology",
        description="File path for the toplogy file for the trajectory",
        subtype="FILE_PATH",
        maxlen=0,
    )
    trajectory: StringProperty(  # type: ignore
        name="Trajectory",
        description="File path for the trajectory file for the trajectory",
        subtype="FILE_PATH",
        maxlen=0,
    )
    name: StringProperty(  # type: ignore
        name="Name",
        description="Name of the molecule on import",
        default="NewTrajectory",
        maxlen=0,
    )
    style: EnumProperty(  # type: ignore
        name="Style",
        description="Default style for importing",
        items=STYLE_ITEMS,
        default="spheres",
    )
    setup_nodes: BoolProperty(  # type: ignore
        name="Setup Nodes",
        description="Add nodes to the scene to load the trajectory",
        default=True,
    )

    def execute(self, context):
        traj = trajectory.load(
            top=self.topology,
            traj=self.trajectory,
            name=self.name,
            style=self.style if self.setup_nodes else None,
        )

        context.view_layer.objects.active = traj.object
        context.scene.frame_start = 0
        context.scene.frame_end = int(traj.universe.trajectory.n_frames - 1)

        self.report(
            {"INFO"},
            message=f"Imported '{self.topology}' as {traj.name} "
            f"with {str(traj.universe.trajectory.n_frames)} "
            f"frames from '{self.trajectory}'.",
        )

        return {"FINISHED"}


class MN_OT_Import_OxDNA_Trajectory(TrajectoryImportOperator):
    """
    Blender operator for importing oxDNA trajectories.
    """

    bl_idname = "mn.import_oxdna"

    def execute(self, context):
        trajectory.load_oxdna(top=self.topology, traj=self.trajectory, name=self.name)
        return {"FINISHED"}


CLASSES = [
    MN_OT_Import_Fetch,
    MN_OT_Import_OxDNA_Trajectory,
    MN_OT_Import_Trajectory,
    MN_OT_Reload_Trajectory,
    MN_OT_Import_Map,
    MN_OT_Import_Star_File,
    MN_OT_Import_Cell_Pack,
    MN_OT_Import_Protein_Local,
    MN_OT_Import_Molecule,
    MN_FH_Import_Molecule,
]
