from pathlib import Path
import bpy
import databpy
import MDAnalysis as mda
from bpy.props import (
    BoolProperty,
    CollectionProperty,
    EnumProperty,
    IntProperty,
    StringProperty,
)
from bpy.types import Context, Operator
from .. import entities
from ..blender.utils import path_resolve
from ..download import CACHE_DIR, FileDownloadPDBError
from ..entities import Molecule, density, ensemble, trajectory
from ..nodes import nodes
from . import node_info
from .style import STYLE_ITEMS


def _add_node(node_name, context, show_options=False, material="default"):
    """
    Add a node group to the node tree and set the values.

    intended to be called upon button press in the node tree, and not for use in general scripting
    """

    # actually invoke the operator to add a node to the current node tree
    # use_transform=True ensures it appears where the user's mouse is and is currently
    # being moved so the user can place it where they wish
    bpy.ops.node.add_node(
        "INVOKE_DEFAULT", type="GeometryNodeGroup", use_transform=True
    )
    node = context.active_node
    node.node_tree = bpy.data.node_groups[node_name]
    node.width = nodes.NODE_WIDTH
    node.show_options = show_options
    node.label = node_name
    node.name = node_name

    # if added node has a 'Material' input, set it to the default MN material
    nodes.assign_material(node, new_material=material)


class MN_OT_Add_Custom_Node_Group(Operator):
    bl_idname = "mn.add_custom_node_group"
    bl_label = "Add Custom Node Group"
    # bl_description = "Add Molecular Nodes custom node group."
    bl_options = {"REGISTER", "UNDO"}
    node_name: StringProperty(  # type: ignore
        name="node_name", description="", default="", subtype="NONE", maxlen=0
    )
    node_label: StringProperty(name="node_label", default="")  # type: ignore
    node_description: StringProperty(  # type: ignore
        name="node_description",
        description="",
        default="Add MolecularNodes custom node group.",
        subtype="NONE",
    )
    node_link: BoolProperty(name="node_link", default=True)  # type: ignore

    @classmethod
    def description(cls, context, properties):
        return properties.node_description

    def execute(self, context):
        try:
            nodes.append(self.node_name, link=self.node_link)
            _add_node(self.node_name, context)  # , label=self.node_label)
        except RuntimeError:
            self.report(
                {"ERROR"},
                message="Failed to add node. Ensure you are not in edit mode.",
            )
            return {"CANCELLED"}
        return {"FINISHED"}


class MN_OT_Assembly_Bio(Operator):
    bl_idname = "mn.assembly_bio"
    bl_label = "Build Biological Assembly"
    bl_description = "Adds node to build biological assembly based on symmetry operations that are extraced from the structure file"
    bl_options = {"REGISTER", "UNDO"}

    inset_node: BoolProperty(default=False)  # type: ignore

    @classmethod
    def poll(cls, context):
        # this just checks to see that there is some biological assembly information that
        # is associated with the object / molecule. If there isn't then the assembly
        # operator will be greyed out and unable to be executed
        obj = context.active_object
        if obj is None:
            return False
        return obj.mn.biological_assemblies != ""

    def execute(self, context):
        obj = context.active_object
        if not isinstance(obj, bpy.types.Object):
            self.report({"ERROR"}, "No active object")
            return {"CANCELLED"}

        with databpy.nodes.DuplicatePrevention():
            try:
                if self.inset_node:
                    nodes.assembly_insert(obj)
                else:
                    tree_assembly = nodes.assembly_initialise(obj)
                    _add_node(tree_assembly.name, context)
            except (KeyError, ValueError) as e:
                self.report({"ERROR"}, "Unable to build biological assembly node.")
                self.report({"ERROR"}, str(e))
                return {"CANCELLED"}

        return {"FINISHED"}


class MN_OT_iswitch_custom(Operator):
    bl_idname = "mn.iswitch_custom"
    # bl_idname = "mn.selection_custom"
    bl_label = "Custom ISwitch Node"
    bl_options = {"REGISTER", "UNDO"}

    description: StringProperty(name="Description")  # type: ignore
    dtype: EnumProperty(  # type: ignore
        name="Data type",
        items=(
            ("RGBA", "RGBA", "Color iswitch."),
            ("BOOLEAN", "BOOLEAN", "Boolean iswitch"),
        ),
    )
    field: StringProperty(name="field", default="chain_id")  # type: ignore
    prefix: StringProperty(name="prefix", default="Chain ")  # type: ignore
    node_property: StringProperty(name="node_property", default="chain_ids")  # type: ignore
    node_name: StringProperty(name="node_name", default="chain")  # type: ignore
    starting_value: IntProperty(name="starting_value", default=0)  # type: ignore

    @classmethod
    def poll(cls, context: Context) -> bool:
        obj = context.active_object
        if obj is None:
            return False
        return True

    @classmethod
    def description(cls, context, properties):
        return properties.description

    def execute(self, context):
        object = context.view_layer.objects.active
        prop = object[self.node_property]
        name = object.name
        if not prop:
            self.report(
                {"WARNING"},
                message=f"{self.node_property} not available for {object.name}.",
            )
            return {"CANCELLED"}

        prefix = {"BOOLEAN": "Select", "RGBA": "Color"}[self.dtype]
        node_name = " ".join([prefix, self.node_name, name])

        with databpy.nodes.DuplicatePrevention():
            node_chains = nodes.custom_iswitch(
                name=node_name,
                dtype=self.dtype,
                iter_list=prop,
                start=self.starting_value,
                field=self.field,
                prefix=self.prefix,
            )

        _add_node(node_chains.name, context)

        return {"FINISHED"}


class MN_OT_Residues_Selection_Custom(Operator):
    bl_idname = "mn.residues_selection_custom"
    bl_label = "Res ID Custom"
    bl_description = "Create a selection based on the provided residue strings.\nThis \
        node is built on a per-molecule basis, taking into account the residues that \
        were input."
    bl_options = {"REGISTER", "UNDO"}

    input_resid_string: StringProperty(  # type: ignore
        name="Select residue IDs: ",
        description="Enter a string value.",
        default="19,94,1-16",
    )

    def execute(self, context):
        with databpy.nodes.DuplicatePrevention():
            node_residues = nodes.resid_multiple_selection(
                node_name="MN_select_res_id_custom",
                input_resid_string=self.input_resid_string,
            )

        _add_node(node_residues.name, context)
        return {"FINISHED"}

    def invoke(self, context, event):
        return context.window_manager.invoke_props_dialog(self)


def get_swap_items(self, context):
    node = context.active_node
    prefix = node.node_tree.name.split(" ")[0].lower()
    if prefix == "is":
        prefix = "select"

    items = [
        (item.name, item.label, item.description)
        for item in node_info.menu_items.get_submenu(prefix).items
        if (not item.is_break and not item.is_custom and item.name != "Set Color")
    ]
    return items


class MN_OT_Node_Swap(Operator):
    bl_idname = "mn.node_swap"
    bl_label = "Swap Node"
    bl_description = "Swap this node for another."

    node_description: StringProperty(default="Swap selected node for another")  # type: ignore
    node_items: EnumProperty(items=get_swap_items)  # type: ignore

    @classmethod
    def description(cls, context, properties):
        return properties.node_description

    def execute(self, context: Context):
        node = context.active_node
        nodes.swap(node, self.node_items)
        return {"FINISHED"}


class MN_OT_Change_Color(Operator):
    bl_idname = "mn.change_color"
    bl_label = "Color"

    color: EnumProperty(  # type: ignore
        items=(
            (item.name, item.label, item.description)
            for item in node_info.menu_items.get_submenu("color").items
            if (not item.is_break and not item.is_custom and item.name != "Set Color")
        )
    )

    def execute(self, context: Context):
        node = context.active_node
        nodes.swap(node, self.color)
        self.report({"INFO"}, f"Selected {self.color}")
        return {"FINISHED"}


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
    remove_solvent: BoolProperty(  # type: ignore
        default=True,
        name="Remove Solvent",
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
        layout.prop(self, "remove_solvent")
        layout.prop(self, "assembly")

        return layout


class MN_OT_Import_Molecule(Import_Molecule):
    bl_idname = "mn.import_molecule"
    bl_label = "Import a Molecule"
    bl_options = {"REGISTER", "UNDO"}

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
                Molecule.load(
                    Path(self.directory, file.name),
                    name=file.name,
                    remove_solvent=self.remove_solvent,
                ).add_style(style, assembly=self.assembly)
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
    bl_description = "Download and open a structure from the PDB"
    bl_options = {"REGISTER", "UNDO"}

    code: StringProperty(  # type: ignore
        name="PDB",
        description="The 4-character PDB code to download",
        options={"TEXTEDIT_UPDATE"},
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
        default=str(CACHE_DIR),
        subtype="DIR_PATH",
    )
    remove_solvent: BoolProperty(  # type: ignore
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

    def execute(self, context):
        try:
            mol = (
                entities.Molecule.fetch(
                    code=self.code,
                    cache=self.cache_dir,
                    format=self.file_format,
                    remove_solvent=self.remove_solvent,
                    database=self.database,
                )
                .add_style(
                    style=self.style if self.node_setup else None,  # type: ignore
                    assembly=self.assembly,
                )
                .centre_molecule(self.centre_type if self.centre else None)
            )

        except FileDownloadPDBError as e:
            self.report({"ERROR"}, str(e))
            if self.file_format == "pdb":
                self.report(
                    {"ERROR"},
                    "There may not be a `.pdb` formatted file available - try a different download format.",
                )
            return {"CANCELLED"}

        message = f"Downloaded {self.code} as {mol.name}"
        try:
            bpy.context.view_layer.objects.active = mol.object  # type: ignore
        except RuntimeError:
            message += " - MolecularNodes collection is disabled"

        self.report({"INFO"}, message=message)

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
        mol = (
            Molecule.load(
                file_path=path_resolve(self.filepath),
                remove_solvent=self.remove_solvent,
            )
            .centre_molecule(self.centre_type if self.centre else None)
            .add_style(
                style=self.style if self.node_setup else None,  # type: ignore
                assembly=self.assembly,
            )
        )

        message = f"Imported '{self.filepath}' as {mol.name}"
        try:
            bpy.context.view_layer.objects.active = mol.object  # type: ignore
        except RuntimeError:
            message += " - MolecularNodes collection is disabled"

        self.report({"INFO"}, message=message)
        return {"FINISHED"}


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
            file_path=path_resolve(self.filepath),
            node_setup=self.node_setup,
        )
        return {"FINISHED"}


class MN_OT_Import_Cell_Pack(ImportEnsemble):
    bl_idname = "mn.import_cell_pack"
    bl_label = "Load"
    bl_description = "Load a CellPack ensemble from a .cif or .bcif file"
    bl_options = {"REGISTER"}

    def execute(self, context):
        ensemble.load_cellpack(
            file_path=path_resolve(self.filepath),
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
            file_path=path_resolve(scene.mn.import_density),
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
        topo = path_resolve(obj.mn.filepath_topology)
        traj = path_resolve(obj.mn.filepath_trajectory)

        if "oxdna" in obj.mn.entity_type:
            uni = mda.Universe(
                topo,
                traj,
                topology_format=trajectory.oxdna.OXDNAParser,
                format=trajectory.oxdna.OXDNAReader,
            )
            traj = trajectory.oxdna.OXDNA(uni)
        else:
            u = mda.Universe(topo, traj)
            traj = trajectory.Trajectory(u)

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
            top=path_resolve(self.topology),
            traj=path_resolve(self.trajectory),
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
        trajectory.load_oxdna(
            top=path_resolve(self.topology),
            traj=path_resolve(self.trajectory),
            name=self.name,
        )
        return {"FINISHED"}


CLASSES = [
    MN_OT_Add_Custom_Node_Group,
    MN_OT_Residues_Selection_Custom,
    MN_OT_Assembly_Bio,
    MN_OT_iswitch_custom,
    MN_OT_Change_Color,
    MN_OT_Node_Swap,
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
