import bpy
from .wwpdb import fetch
from .retrieve import FileDownloadPDBError

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


class MN_OT_Import_AlphaFold(bpy.types.Operator):
    bl_idname = "mn.import_alphafold"
    bl_label = "Fetch"
    bl_description = "Download specified structure from the AlphaFold databse"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        scene = context.scene
        pdb_code = scene.MN_alphafold_code
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


def panel(layout, scene):
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
