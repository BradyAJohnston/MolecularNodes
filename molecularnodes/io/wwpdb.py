import bpy
from pathlib import Path
from . import parse
from .retrieve import download


def fetch(
    pdb_code,
    style='spheres',
    centre=False,
    del_solvent=True,
    cache_dir=None,
    build_assembly=False,
    format="mmtf"
):

    if build_assembly:
        centre = False

    file_path = download(code=pdb_code, format=format, cache=cache_dir)

    parsers = {
        'mmtf': parse.MMTF,
        'pdb': parse.PDB,
        'cif': parse.CIF
    }
    molecule = parsers[format](file_path=file_path)

    model = molecule.create_model(
        name=pdb_code,
        centre=centre,
        style=style,
        del_solvent=del_solvent,
        build_assembly=build_assembly
    )

    model.mn['pdb_code'] = pdb_code
    model.mn['molecule_type'] = 'pdb'

    return molecule

# Properties that can be set in the scene, to be passed to the operator


bpy.types.Scene.MN_pdb_code = bpy.props.StringProperty(
    name='PDB',
    description='The 4-character PDB code to download',
    options={'TEXTEDIT_UPDATE'},
    maxlen=4
)
bpy.types.Scene.MN_cache_dir = bpy.props.StringProperty(
    name='',
    description='Directory to save the downloaded files',
    options={'TEXTEDIT_UPDATE'},
    default=str(Path('~', '.MolecularNodes').expanduser()),
    subtype='DIR_PATH'
)
bpy.types.Scene.MN_cache = bpy.props.BoolProperty(
    name="Cache Downloads",
    description="Save the downloaded file in the given directory",
    default=True
)
bpy.types.Scene.MN_import_format_download = bpy.props.EnumProperty(
    name="Format",
    description="Format to download as from the PDB",
    items=(
        ("mmtf", ".mmtf", "The binary compressed MMTF, fastest for downloading"),
        ("cif", ".cif", 'The new standard of .cif / .mmcif'),
        ("pdb", ".pdb", "The classic (and depcrecated) PDB format")
    )
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

        if not scene.MN_cache:
            cache_dir = None

        style = None
        if scene.MN_import_node_setup:
            style = scene.MN_import_style

        mol = fetch(
            pdb_code=pdb_code,
            centre=scene.MN_import_centre,
            del_solvent=scene.MN_import_del_solvent,
            style=style,
            cache_dir=cache_dir,
            build_assembly=scene.MN_import_build_assembly,
            format=scene.MN_import_format_download
        )

        bpy.context.view_layer.objects.active = mol.object
        self.report(
            {'INFO'}, message=f"Imported '{pdb_code}' as {mol.object.name}")

        return {"FINISHED"}

# the UI for the panel, which will display the operator and the properties


def panel(layout, scene):

    layout.label(text="Download from PDB", icon="IMPORT")
    layout.separator()
    row_import = layout.row().split(factor=0.5)
    row_import.prop(scene, 'MN_pdb_code')
    download = row_import.split(factor=0.3)
    download.prop(scene, 'MN_import_format_download', text="")
    download.operator('mn.import_wwpdb')
    layout.separator(factor=0.4)
    row = layout.row().split(factor=0.3)
    row.prop(scene, 'MN_cache')
    row_cache = row.row()
    row_cache.prop(scene, 'MN_cache_dir')
    row_cache.enabled = scene.MN_cache
    layout.separator()
    layout.label(text="Options", icon="MODIFIER")
    options = layout.column(align=True)
    row = options.row()
    row.prop(scene, 'MN_import_node_setup', text="")
    col = row.column()
    col.prop(scene, "MN_import_style")
    col.enabled = scene.MN_import_node_setup

    options.separator()
    grid = options.grid_flow()
    grid.prop(scene, 'MN_import_build_assembly')
    grid.prop(scene, 'MN_import_centre')
    grid.prop(scene, 'MN_import_del_solvent')
