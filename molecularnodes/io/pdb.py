import bpy
from pathlib import Path
from . import parse

bpy.types.Scene.MN_pdb_code = bpy.props.StringProperty(
    name = 'PDB', 
    description = 'The 4-character PDB code to download', 
    options = {'TEXTEDIT_UPDATE'}, 
    maxlen = 4
    )
bpy.types.Scene.MN_cache_dir = bpy.props.StringProperty(
    name = 'Cache',
    description = 'Location to cache PDB files',
    options = {'TEXTEDIT_UPDATE'},
    default = str(Path('~', '.MolecularNodes').expanduser()),
    subtype = 'FILE_PATH'
)
bpy.types.Scene.MN_cache = bpy.props.BoolProperty(
    name = "Cache", 
    default = True
)
bpy.types.Scene.MN_import_format_download = bpy.props.EnumProperty(
    name = "Format",
    description = "Format to download as from the PDB",
    items = (
        ("mmtf", ".mmtf", "The binary compressed MMTF, fastest for downloading."),
        ("pdb", ".pdb", "The classic (and depcrecated) PDB format"), 
        ("cif", ".mmcif", 'The new standard of .mmcif')
    )
)

def load(
    pdb_code,             
    style = 'spheres',               
    centre = False,               
    del_solvent = True,
    cache_dir = None,
    build_assembly = False
    ):
    import biotite.database.rcsb as rcsb
    
    if build_assembly:
        centre = False
        
    file_format = "mmtf"
    file_path = rcsb.fetch(pdb_code, file_format, target_path=cache_dir)
    
    datafile = parse.MMTF(file_path=file_path)
    
    mol = datafile.create_model(
        name=pdb_code, 
        centre=centre, 
        del_solvent=del_solvent
    )
    
    mol.mn['pdb_code'] = pdb_code
    mol.mn['molecule_type'] = 'pdb'
    
    return mol

# operator that calls the function to import the structure from the PDB
class MN_OT_Import_Protein_RCSB(bpy.types.Operator):
    bl_idname = "mn.import_protein_rcsb"
    bl_label = "import_protein_fetch_pdb"
    bl_description = "Download and open a structure from the Protein Data Bank"
    bl_options = {"REGISTER", "UNDO"}

    @classmethod
    def poll(cls, context):
        return not False

    def execute(self, context):
        scene = context.scene
        pdb_code = scene.MN_pdb_code
        cache_dir = scene.MN_cache_dir
        
        if not scene.MN_cache:
            cache_dir = None
        
        style = None
        if scene.MN_import_node_setup:
            style = scene.MN_import_style
        
        mol = load(
            pdb_code=pdb_code,
            centre=scene.MN_import_centre, 
            del_solvent=scene.MN_import_del_solvent,
            style=style,
            cache_dir=cache_dir, 
            node_setup=scene.MN_import_node_setup,
            build_assembly = scene.MN_import_build_assembly
        )
        
        bpy.context.view_layer.objects.active = mol
        self.report({'INFO'}, message=f"Imported '{pdb_code}' as {mol.name}")
        
        return {"FINISHED"}

def panel(layout, scene):
    
    layout.label(text = "Download from PDB", icon="IMPORT")
    layout.separator()
    row_import = layout.row()
    row_import.prop(scene, 'MN_import_format_download')
    row_import.prop(scene, 'MN_pdb_code')
    row_import.operator('mn.import_protein_rcsb', text='Download')
    layout.separator()
    layout.label(text = "Options", icon = "MODIFIER")
    options = layout.column(align = True)
    row = options.row()
    row.prop(scene, 'MN_import_node_setup', text = "")
    col = row.column()
    col.prop(scene, "MN_import_style")
    col.enabled = scene.MN_import_node_setup
    
    options.separator()
    row = options.row()
    row.prop(scene, 'MN_cache', text="")
    col = row.column()
    col.prop(scene, 'MN_cache_dir', text = "Cache")
    col.enabled = scene.MN_cache
    grid = options.grid_flow()
    grid.prop(scene, 'MN_import_build_assembly')
    grid.prop(scene, 'MN_import_centre')
    grid.prop(scene, 'MN_import_del_solvent')