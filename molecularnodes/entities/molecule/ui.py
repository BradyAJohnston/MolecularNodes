from pathlib import Path

import bpy
from bpy.types import Context, UILayout
from bpy.props import EnumProperty, StringProperty, BoolProperty, CollectionProperty
from biotite import InvalidFileError
import os
import io

from ...download import FileDownloadPDBError, download, CACHE_DIR
from ...blender import path_resolve
from .oldcif import OldCIF
from .molecule import Molecule
from .pdb import PDB
from .pdbx import BCIF, CIF
from .sdf import SDF
from ...style import STYLE_ITEMS


def parse(filepath) -> Molecule:
    # TODO: I don't like that we might be dealing with bytes or a filepath here,
    # I need to work out a nicer way to have it be cleanly one or the other

    if isinstance(filepath, io.BytesIO):
        suffix = ".bcif"
    else:
        filepath = path_resolve(filepath)
        suffix = Path(filepath).suffix

    parser = {
        ".pdb": PDB,
        ".pdbx": CIF,
        ".cif": CIF,
        ".bcif": BCIF,
        ".mol": SDF,
        ".sdf": SDF,
    }

    if suffix not in parser:
        raise ValueError(f"Unable to open local file. Format '{suffix}' not supported.")
    try:
        molecule = parser[suffix](filepath)
    except InvalidFileError:
        molecule = OldCIF(filepath)

    return molecule


def fetch(
    pdb_code: str,
    style: str | None = "spheres",
    centre: str = "",
    del_solvent: bool = True,
    del_hydrogen: bool = False,
    cache_dir: str | None = None,
    build_assembly: bool = False,
    database: str = "rcsb",
    format: str = "bcif",
    color: str = "common",
) -> Molecule:
    if build_assembly:
        centre = ""

    file_path = download(
        code=pdb_code, format=format, cache=cache_dir, database=database
    )

    mol = parse(file_path)

    obj = mol.create_object(
        name=pdb_code,
        centre=centre,
        style=style,
        del_solvent=del_solvent,
        del_hydrogen=del_hydrogen,
        build_assembly=build_assembly,
        color=color,
    )

    obj.mn["pdb_code"] = pdb_code
    obj.mn["entity_type"] = format

    return mol


def load_local(
    file_path,
    name="Name",
    centre="",
    del_solvent=True,
    del_hydrogen=False,
    style="spheres",
    build_assembly=False,
):
    mol = parse(file_path)
    mol.create_object(
        name=name,
        style=style,
        build_assembly=build_assembly,
        centre=centre,
        del_solvent=del_solvent,
        del_hydrogen=del_hydrogen,
    )
    return mol


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
    centre: EnumProperty(  # type: ignore
        name="Centre",
        description="Centre the structure at the world origin using the given method",
        default="None",
        items=(
            ("None", "None", "No centering is applied", 1),
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

    def draw(self, context: Context) -> UILayout:
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

    def draw(self, context: Context):
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


class MN_OT_Import_Fetch(Import_Molecule):
    bl_idname = "mn.import_fetch"
    bl_label = "Download a Molecule"
    bl_description = "Download a molecule from the wwPDB and import it to the scene"

    code: StringProperty(  # type: ignore
        default="4ozs",
        name="PDB Code",
        description="Code to use for downloading from the wwPDB",
    )
    file_format: EnumProperty(  # type: ignore
        name="Format",
        description="Format to download as from the PDB",
        default="bcif",
        items=DOWNLOAD_FORMATS,
    )

    database: EnumProperty(  # type: ignore
        name="Database",
        description="Database to download structure from",
        default="rcsb",
        items=(
            ("rcsb", "rcsb", "The RCSB PDB"),
            ("pdbe", "PDBe", "The european protein data bank"),
            ("pdbj", "PDBj", "The Japanese protein data bank"),
            ("alphafold", "AplhaFold", "The AlphaFold database"),
            ("emdb", "EM Database", "Electron microscopy database"),
        ),
    )

    def execute(self, context):
        try:
            file_path = download(
                self.code, format=self.file_format, cache=self.cache, database="rcsb"
            )
            mol = parse(file_path)
            mol.create_object(
                name=self.code,
                style=self.style,
                centre=self.centre,
                del_solvent=self.del_solvent,
                build_assembly=self.assembly,
            )
        except Exception:
            return {"CANCELLED"}
        return {"FINISHED"}

    def invoke(self, context, event):
        context.window_manager.invoke_props_dialog(self)
        return {"RUNNING_MODAL"}


# Properties that can be set in the scene, to be passed to the operator


bpy.types.Scene.MN_pdb_code = StringProperty(
    name="PDB",
    description="The 4-character PDB code to download",
    options={"TEXTEDIT_UPDATE"},
    maxlen=4,
)
bpy.types.Scene.MN_cache_dir = StringProperty(
    name="",
    description="Directory to save the downloaded files",
    options={"TEXTEDIT_UPDATE"},
    default=CACHE_DIR,
    subtype="DIR_PATH",
)
bpy.types.Scene.MN_cache = BoolProperty(
    name="Cache Downloads",
    description="Save the downloaded file in the given directory",
    default=True,
)
bpy.types.Scene.MN_import_del_hydrogen = BoolProperty(
    name="Remove Hydrogens",
    description="Remove the hydrogens from a structure on import",
    default=False,
)
bpy.types.Scene.MN_import_format_download = EnumProperty(
    name="Format",
    description="Format to download as from the PDB",
    items=(
        ("bcif", ".bcif", "Binary compressed .cif file, fastest for downloading"),
        ("cif", ".cif", "The new standard of .cif / .mmcif"),
        ("pdb", ".pdb", "The classic (and depcrecated) PDB format"),
    ),
)
bpy.types.Scene.MN_import_local_path = StringProperty(
    name="File",
    description="File path of the structure to open",
    options={"TEXTEDIT_UPDATE"},
    subtype="FILE_PATH",
    maxlen=0,
)
bpy.types.Scene.MN_import_local_name = StringProperty(
    name="Name",
    description="Name of the molecule on import",
    options={"TEXTEDIT_UPDATE"},
    default="NewMolecule",
    maxlen=0,
)

bpy.types.Scene.MN_alphafold_code = StringProperty(
    name="UniProt ID",
    description="The UniProt ID to use for downloading from the AlphaFold databse",
    options={"TEXTEDIT_UPDATE"},
)

bpy.types.Scene.MN_import_format_alphafold = EnumProperty(
    name="Format",
    description="Format to download as from the PDB",
    items=(
        # ("bcif", ".bcif", "Binary compressed .cif file, fastest for downloading"),
        ("cif", ".cif", "The new standard of .cif / .mmcif"),
        ("pdb", ".pdb", "The classic (and depcrecated) PDB format"),
    ),
)


# operator that is called by the 'button' press which calls the fetch function


class MN_OT_Import_wwPDB(bpy.types.Operator):
    bl_idname = "mn.import_wwpdb"
    bl_label = "Fetch"
    bl_description = "Download and open a structure from the Protein Data Bank"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        scene = context.scene
        pdb_code = scene.MN_pdb_code
        cache_dir = scene.MN_cache_dir
        file_format = scene.MN_import_format_download

        if not scene.MN_cache:
            cache_dir = None

        style = None
        if scene.mn.import_node_setup:
            style = scene.mn.import_style

        centre = ""
        if scene.mn.import_centre:
            centre = scene.mn.centre_type

        try:
            mol = fetch(
                pdb_code=pdb_code,
                centre=centre,
                del_solvent=scene.mn.import_del_solvent,
                del_hydrogen=scene.MN_import_del_hydrogen,
                style=style,
                cache_dir=cache_dir,
                build_assembly=scene.mn.import_build_assembly,
                format=file_format,
            )
        except FileDownloadPDBError as e:
            self.report({"ERROR"}, str(e))
            if file_format == "pdb":
                self.report(
                    {"ERROR"},
                    "There may not be a `.pdb` formatted file available - try a different download format.",
                )
            return {"CANCELLED"}

        bpy.context.view_layer.objects.active = mol.object
        self.report({"INFO"}, message=f"Imported '{pdb_code}' as {mol.object.name}")

        return {"FINISHED"}


class MN_OT_Import_Protein_Local(bpy.types.Operator):
    bl_idname = "mn.import_protein_local"
    bl_label = "Load"
    bl_description = "Open a local structure file"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        scene = context.scene
        file_path = scene.MN_import_local_path

        style = scene.mn.import_style
        if not scene.mn.import_node_setup:
            style = None

        centre = ""
        if scene.mn.import_centre:
            centre = scene.mn.centre_type

        mol = load_local(
            file_path=file_path,
            name=scene.MN_import_local_name,
            centre=centre,
            del_solvent=scene.mn.import_del_solvent,
            del_hydrogen=scene.MN_import_del_hydrogen,
            style=style,
            build_assembly=scene.mn.import_build_assembly,
        )

        # return the good news!
        bpy.context.view_layer.objects.active = mol.object
        self.report({"INFO"}, message=f"Imported '{file_path}' as {mol.name}")
        return {"FINISHED"}

    def invoke(self, context, event):
        return self.execute(context)


class MN_OT_Import_AlphaFold(bpy.types.Operator):
    bl_idname = "mn.import_alphafold"
    bl_label = "Fetch"
    bl_description = "Download specified structure from the AlphaFold databse"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        scene = context.scene
        pdb_code = scene.MN_alphafold_code.strip()
        cache_dir = scene.MN_cache_dir
        file_format = scene.MN_import_format_alphafold

        if not scene.MN_cache:
            cache_dir = None

        style = None
        if scene.mn.import_node_setup:
            style = scene.mn.import_style

        centre = ""
        if scene.mn.import_centre:
            centre = scene.mn.centre_type

        try:
            mol = fetch(
                pdb_code=pdb_code,
                centre=centre,
                del_solvent=scene.mn.import_del_solvent,
                style=style,
                cache_dir=cache_dir,
                build_assembly=scene.mn.import_build_assembly,
                format=file_format,
                database="alphafold",
                color="plddt",
            )
        except FileDownloadPDBError as e:
            self.report({"ERROR"}, str(e))
            if file_format == "pdb":
                self.report(
                    {"ERROR"},
                    "There may not be a `.pdb` formatted file available - try a different download format.",
                )
            return {"CANCELLED"}

        bpy.context.view_layer.objects.active = mol.object
        self.report({"INFO"}, message=f"Imported '{pdb_code}' as {mol.object.name}")

        return {"FINISHED"}


# the UI for the panel, which will display the operator and the properties


def check_online_access_for_ui(layout: bpy.types.UILayout) -> bpy.types.UILayout:
    """
    Disable UI without Online Access

    Checks for the online access permissions, and adds a warning and disables following
    UI elements if it fails the check. Returns the UILayout that will have .enabled flag
    set to False, disabling all subsequent uses of the layout.

    Args:
        layout (bpy.types.UILayout): The UILayout element to add the warning and potentially
        disable.

    Returns:
        bpy.types.UILayout: The altered UILayout element, for use in downstream UI
        components.
    """
    if not bpy.app.online_access:
        layout.label(
            text="Online access disabled. Change in Blender's system preferences.",
            icon="ERROR",
        )
        op = layout.operator("wm.url_open", text="Online Access Docs", icon="URL")
        op.url = "https://docs.blender.org/manual/en/dev/editors/preferences/system.html#bpy-types-preferencessystem-use-online-access"
        layout = layout.column()
        layout.alert = True
        layout.enabled = False

    return layout


def panel_wwpdb(layout, scene):
    layout.label(text="Download from PDB", icon="IMPORT")
    layout.separator()

    layout = check_online_access_for_ui(layout)

    row_import = layout.row().split(factor=0.5)
    row_import.prop(scene, "MN_pdb_code")
    download = row_import.split(factor=0.3)
    download.prop(scene, "MN_import_format_download", text="")
    download.operator("mn.import_wwpdb")
    layout.separator(factor=0.4)

    row = layout.row().split(factor=0.3)
    row.prop(scene, "MN_cache")
    row_cache = row.row()
    row_cache.prop(scene, "MN_cache_dir")
    row_cache.enabled = scene.MN_cache
    layout.separator()

    layout.label(text="Options", icon="MODIFIER")
    options = layout.column(align=True)

    row = options.row()
    row.prop(scene.mn, "import_node_setup", text="")
    col = row.column()
    col.prop(scene.mn, "import_style")
    col.enabled = scene.mn.import_node_setup

    row_centre = options.row()
    row_centre.prop(scene.mn, "import_centre", icon_value=0)
    col_centre = row_centre.column()
    col_centre.prop(scene.mn, "centre_type", text="")
    col_centre.enabled = scene.mn.import_centre
    options.separator()

    grid = options.grid_flow()
    grid.prop(scene.mn, "import_build_assembly")
    grid.prop(scene.mn, "import_del_solvent")
    grid.prop(scene, "MN_import_del_hydrogen")


def panel_alphafold(layout, scene):
    layout.label(text="Download from the AlphaFold DataBase", icon="IMPORT")
    layout.separator()

    layout = check_online_access_for_ui(layout)

    row_import = layout.row().split(factor=0.5)
    row_import.prop(scene, "MN_alphafold_code")
    download = row_import.split(factor=0.3)
    download.prop(scene, "MN_import_format_alphafold", text="")
    download.operator("mn.import_alphafold")
    layout.separator(factor=0.4)

    row = layout.row().split(factor=0.3)
    row.prop(scene, "MN_cache")
    row_cache = row.row()
    row_cache.prop(scene, "MN_cache_dir")
    row_cache.enabled = scene.MN_cache
    layout.separator()

    layout.label(text="Options", icon="MODIFIER")
    options = layout.column(align=True)

    row = options.row()
    row.prop(scene.mn, "import_node_setup", text="")
    col = row.column()
    col.prop(scene.mn, "import_style")
    col.enabled = scene.mn.import_node_setup

    row_centre = options.row()
    row_centre.prop(scene.mn, "import_centre", icon_value=0)
    col_centre = row_centre.column()
    col_centre.prop(scene.mn, "centre_type", text="")
    col_centre.enabled = scene.mn.import_centre
    options.separator()

    grid = options.grid_flow()
    grid.prop(scene.mn, "import_build_assembly")
    grid.prop(scene.mn, "import_del_solvent")
    # grid.prop(scene, "MN_import_del_hydrogen")


# operator that calls the function to import the structure from a local file


def panel_local(layout, scene):
    layout.label(text="Load a Local File", icon="FILE_TICK")
    layout.separator()

    row_name = layout.row(align=False)
    row_name.prop(scene, "MN_import_local_name")
    row_name.operator("mn.import_protein_local")

    row_import = layout.row()
    row_import.prop(scene, "MN_import_local_path")
    layout.separator()

    layout.label(text="Options", icon="MODIFIER")
    options = layout.column(align=True)

    row = options.row()
    row.prop(scene.mn, "import_node_setup", text="")
    col = row.column()
    col.prop(scene.mn, "import_style")
    col.enabled = scene.mn.import_node_setup

    row_centre = options.row()

    row_centre.prop(scene.mn, "import_centre", icon_value=0)
    # row_centre.prop()
    col_centre = row_centre.column()
    col_centre.prop(scene.mn, "centre_type", text="")
    col_centre.enabled = scene.mn.import_centre
    options.separator()

    grid = options.grid_flow()
    grid.prop(scene.mn, "import_build_assembly")
    grid.prop(scene.mn, "import_del_solvent", icon_value=0)
    grid.prop(scene, "MN_import_del_hydrogen", icon_value=0)


CLASSES = [
    MN_OT_Import_AlphaFold,
    MN_OT_Import_Protein_Local,
    MN_OT_Import_wwPDB,
    MN_OT_Import_Molecule,
    MN_FH_Import_Molecule,
    MN_OT_Import_Fetch,
]
