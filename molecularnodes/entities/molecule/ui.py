from pathlib import Path

import bpy
from bpy.types import Context, UILayout
from bpy.props import EnumProperty, StringProperty, BoolProperty, CollectionProperty
from bpy_extras.io_utils import ImportHelper
from biotite import InvalidFileError
import os
import io
from pathlib import Path

from ...download import FileDownloadPDBError, download, CACHE_DIR
from ...blender import path_resolve
from .molecule import Molecule
from .pdb import PDB
from .pdbx import BCIF, CIF
from .sdf import SDF
from ...style import STYLE_ITEMS


def addon_preferences() -> bpy.types.AddonPreferences:
    try:
        return bpy.context.preferences.addons[__package__].preferences
    except KeyError:
        return bpy.context.preferences.addons[
            "bl_ext.vscode_development.molecularnodes"
        ].preferences


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

    return parser[suffix](filepath)


def fetch(
    code: str,
    style: str | None = "spheres",
    centre: str | None = None,
    del_solvent: bool = True,
    del_hydrogen: bool = False,
    cache_dir: str | Path | None = None,
    build_assembly: bool = False,
    database: str = "rcsb",
    format: str = "bcif",
    color: str = "common",
) -> Molecule:
    if build_assembly:
        centre = ""

    file_path = download(code=code, format=format, cache=cache_dir, database=database)

    mol = parse(file_path)

    obj = mol.create_object(
        name=code,
        centre=centre,
        style=style,
        del_solvent=del_solvent,
        del_hydrogen=del_hydrogen,
        build_assembly=build_assembly,
        color=color,
    )

    obj.mn["code"] = code
    obj.mn["entity_type"] = format

    return mol


def load_local(
    file_path,
    centre: str | None = "",
    style: str | None = "spheres",
    del_solvent=True,
    del_hydrogen=False,
    build_assembly=False,
):
    mol = parse(file_path)
    mol.create_object(
        name=Path(file_path).stem,
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
            mol = fetch(
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
        mol = parse(self.filepath)
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
    row_import.prop(scene.mn, "import_code_pdb")
    download = row_import.split(factor=0.3)
    download.prop(scene.mn, "import_format_wwpdb", text="")
    op = download.operator("mn.import_fetch")
    op.code = scene.mn.import_code_pdb
    op.database = "wwpdb"
    op.file_format = scene.mn.import_format_wwpdb
    op.node_setup = scene.mn.import_node_setup
    op.assembly = scene.mn.import_build_assembly
    op.style = scene.mn.import_style
    op.centre = scene.mn.import_centre
    op.centre_type = scene.mn.import_centre_type
    op.cache_dir = str(addon_preferences().cache_dir)
    layout.separator(factor=0.4)

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
    col_centre.prop(scene.mn, "import_centre_type", text="")
    col_centre.enabled = scene.mn.import_centre
    options.separator()

    grid = options.grid_flow()
    grid.prop(scene.mn, "import_build_assembly")
    grid.prop(scene.mn, "import_del_solvent")
    grid.prop(scene.mn, "import_del_hydrogen")


def panel_alphafold(layout, scene):
    layout.label(text="Download from the AlphaFold DataBase", icon="IMPORT")
    layout.separator()

    layout = check_online_access_for_ui(layout)

    row_import = layout.row().split(factor=0.5)
    row_import.prop(scene.mn, "import_code_alphafold")
    download = row_import.split(factor=0.3)
    download.prop(scene.mn, "import_format_alphafold", text="")
    op = download.operator("mn.import_fetch")
    op.code = scene.mn.import_code_alphafold
    op.database = "alphafold"
    op.file_format = scene.mn.import_format_alphafold
    op.node_setup = scene.mn.import_node_setup
    op.assembly = scene.mn.import_build_assembly
    op.style = scene.mn.import_style
    op.centre = scene.mn.import_centre
    op.centre_type = scene.mn.import_centre_type
    op.cache_dir = str(addon_preferences().cache_dir)

    layout.separator(factor=0.4)

    row = layout.row().split(factor=0.3)
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


# operator that calls the function to import the structure from a local file


def panel_local(layout, scene):
    layout.label(text="Load a Local File", icon="FILE_TICK")
    layout.separator()

    row = layout.row()
    row.prop(scene.mn, "import_local_path")
    op = row.operator("mn.import_local")
    op.filepath = scene.mn.import_local_path
    op.node_setup = scene.mn.import_node_setup
    op.assembly = scene.mn.import_build_assembly
    op.style = scene.mn.import_style
    op.centre = scene.mn.import_centre
    op.centre_type = scene.mn.import_centre_type
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
    col_centre.prop(scene.mn, "import_centre_type", text="")
    col_centre.enabled = scene.mn.import_centre
    options.separator()

    grid = options.grid_flow()
    grid.prop(scene.mn, "import_build_assembly")
    grid.prop(scene.mn, "import_del_solvent", icon_value=0)
    grid.prop(scene.mn, "import_del_hydrogen", icon_value=0)


CLASSES = [
    MN_OT_Import_Fetch,
    MN_OT_Import_Protein_Local,
    MN_OT_Import_Molecule,
    MN_FH_Import_Molecule,
]
