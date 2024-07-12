from pathlib import Path

import bpy
from biotite import InvalidFileError

from ...download import FileDownloadPDBError, download, CACHE_DIR
from ..ensemble.cif import OldCIF
from .molecule import Molecule
from .pdb import PDB
from .pdbx import BCIF, CIF
from .sdf import SDF


def fetch(
    pdb_code,
    style="spheres",
    centre="",
    del_solvent=True,
    cache_dir=None,
    build_assembly=False,
    database: str = "rcsb",
    format="bcif",
    color="common",
) -> Molecule:
    if build_assembly:
        centre = ""

    file_path = download(
        code=pdb_code, format=format, cache=cache_dir, database=database
    )

    parsers = {"pdb": PDB, "cif": CIF, "bcif": BCIF}
    molecule = parsers[format](file_path=file_path)

    model = molecule.create_model(
        name=pdb_code,
        centre=centre,
        style=style,
        del_solvent=del_solvent,
        build_assembly=build_assembly,
        color=color,
    )

    model.mn["pdb_code"] = pdb_code
    model.mn["molecule_type"] = format

    return molecule


def load_local(
    file_path,
    name="Name",
    centre="",
    del_solvent=True,
    style="spheres",
    build_assembly=False,
):
    suffix = Path(file_path).suffix
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
        molecule = parser[suffix](file_path)
    except InvalidFileError:
        molecule = OldCIF(file_path)

    molecule.create_model(
        name=name,
        style=style,
        build_assembly=build_assembly,
        centre=centre,
        del_solvent=del_solvent,
    )
    return molecule


# Properties that can be set in the scene, to be passed to the operator


bpy.types.Scene.MN_pdb_code = bpy.props.StringProperty(
    name="PDB",
    description="The 4-character PDB code to download",
    options={"TEXTEDIT_UPDATE"},
    maxlen=4,
)
bpy.types.Scene.MN_cache_dir = bpy.props.StringProperty(
    name="",
    description="Directory to save the downloaded files",
    options={"TEXTEDIT_UPDATE"},
    default=CACHE_DIR,
    subtype="DIR_PATH",
)
bpy.types.Scene.MN_cache = bpy.props.BoolProperty(
    name="Cache Downloads",
    description="Save the downloaded file in the given directory",
    default=True,
)
bpy.types.Scene.MN_import_format_download = bpy.props.EnumProperty(
    name="Format",
    description="Format to download as from the PDB",
    items=(
        ("bcif", ".bcif", "Binary compressed .cif file, fastest for downloading"),
        ("cif", ".cif", "The new standard of .cif / .mmcif"),
        ("pdb", ".pdb", "The classic (and depcrecated) PDB format"),
    ),
)
bpy.types.Scene.MN_import_local_path = bpy.props.StringProperty(
    name="File",
    description="File path of the structure to open",
    options={"TEXTEDIT_UPDATE"},
    subtype="FILE_PATH",
    maxlen=0,
)
bpy.types.Scene.MN_import_local_name = bpy.props.StringProperty(
    name="Name",
    description="Name of the molecule on import",
    options={"TEXTEDIT_UPDATE"},
    default="NewMolecule",
    maxlen=0,
)
bpy.types.Scene.MN_alphafold_code = bpy.props.StringProperty(
    name="UniProt ID",
    description="The UniProt ID to use for downloading from the AlphaFold databse",
    options={"TEXTEDIT_UPDATE"},
)

bpy.types.Scene.MN_import_format_alphafold = bpy.props.EnumProperty(
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
        if scene.MN_import_node_setup:
            style = scene.MN_import_style

        centre = ""
        if scene.MN_import_centre:
            centre = scene.MN_centre_type

        try:
            mol = fetch(
                pdb_code=pdb_code,
                centre=centre,
                del_solvent=scene.MN_import_del_solvent,
                style=style,
                cache_dir=cache_dir,
                build_assembly=scene.MN_import_build_assembly,
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

        style = scene.MN_import_style
        if not scene.MN_import_node_setup:
            style = None

        centre = ""
        if scene.MN_import_centre:
            centre = scene.MN_centre_type

        mol = load_local(
            file_path=file_path,
            name=scene.MN_import_local_name,
            centre=centre,
            del_solvent=scene.MN_import_del_solvent,
            style=style,
            build_assembly=scene.MN_import_build_assembly,
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
        if scene.MN_import_node_setup:
            style = scene.MN_import_style

        centre = ""
        if scene.MN_import_centre:
            centre = scene.MN_centre_type

        try:
            mol = fetch(
                pdb_code=pdb_code,
                centre=centre,
                del_solvent=scene.MN_import_del_solvent,
                style=style,
                cache_dir=cache_dir,
                build_assembly=scene.MN_import_build_assembly,
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


def panel_wwpdb(layout, scene):
    layout.label(text="Download from PDB", icon="IMPORT")
    layout.separator()

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
    row.prop(scene, "MN_import_node_setup", text="")
    col = row.column()
    col.prop(scene, "MN_import_style")
    col.enabled = scene.MN_import_node_setup

    row_centre = options.row()
    row_centre.prop(scene, "MN_import_centre", icon_value=0)
    col_centre = row_centre.column()
    col_centre.prop(scene, "MN_centre_type", text="")
    col_centre.enabled = scene.MN_import_centre
    options.separator()

    grid = options.grid_flow()
    grid.prop(scene, "MN_import_build_assembly")
    grid.prop(scene, "MN_import_del_solvent")


def panel_alphafold(layout, scene):
    layout.label(text="Download from the AlphaFold DataBase", icon="IMPORT")
    layout.separator()

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
    row.prop(scene, "MN_import_node_setup", text="")
    col = row.column()
    col.prop(scene, "MN_import_style")
    col.enabled = scene.MN_import_node_setup

    row_centre = options.row()
    row_centre.prop(scene, "MN_import_centre", icon_value=0)
    col_centre = row_centre.column()
    col_centre.prop(scene, "MN_centre_type", text="")
    col_centre.enabled = scene.MN_import_centre
    options.separator()

    grid = options.grid_flow()
    grid.prop(scene, "MN_import_build_assembly")
    grid.prop(scene, "MN_import_del_solvent")


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
    row.prop(scene, "MN_import_node_setup", text="")
    col = row.column()
    col.prop(scene, "MN_import_style")
    col.enabled = scene.MN_import_node_setup

    row_centre = options.row()

    row_centre.prop(scene, "MN_import_centre", icon_value=0)
    # row_centre.prop()
    col_centre = row_centre.column()
    col_centre.prop(scene, "MN_centre_type", text="")
    col_centre.enabled = scene.MN_import_centre
    options.separator()

    grid = options.grid_flow()
    grid.prop(scene, "MN_import_build_assembly")
    grid.prop(scene, "MN_import_del_solvent", icon_value=0)


CLASSES = [MN_OT_Import_AlphaFold, MN_OT_Import_Protein_Local, MN_OT_Import_wwPDB]
