import bpy
from pathlib import Path
import numpy as np
from .load import create_molecule
from ..blender import nodes
from .. import assembly

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

def load(
    pdb_code,             
    style = 'spheres',               
    centre = False,               
    del_solvent = True,               
    setup_nodes = True,
    cache_dir = None,
    build_assembly = False
    ):
    from biotite import InvalidFileError
    mol, file = open_structure_rcsb(
        pdb_code = pdb_code, 
        cache_dir = cache_dir
        )
    
    if build_assembly:
        centre = False
    
    mol, coll_frames = create_molecule(
        array = mol,
        name = pdb_code,
        file = file,
        calculate_ss = False,
        centre = centre,
        del_solvent = del_solvent, 
        )
    
    mol.mn['pdb_code'] = pdb_code
    mol.mn['molecule_type'] = 'pdb'
    
    if setup_nodes:
        nodes.create_starting_node_tree(
            object = mol, 
            coll_frames=coll_frames, 
            style = style
            )
    
    try:
        parsed_assembly_file = assembly.mmtf.MMTFAssemblyParser(file)
        mol['biological_assemblies'] = parsed_assembly_file.get_assemblies()
    except InvalidFileError:
        pass
    
    if build_assembly:
        nodes.assembly_insert(mol)
    
    return mol


def get_chain_entity_id(file):
    entities = file['entityList']
    chain_names = file['chainNameList']    
    ent_dic = {}
    for i, ent in enumerate(entities):
        for chain_idx in ent['chainIndexList']:
            chain_id = chain_names[chain_idx]
            if  chain_id in ent_dic.keys():
                next
            else:
                ent_dic[chain_id] = i
    
    return ent_dic

def set_atom_entity_id(mol, file):
    mol.add_annotation('entity_id', int)
    ent_dic = get_chain_entity_id(file)
    
    entity_ids = np.array([ent_dic[x] for x in mol.chain_id])
    
    # entity_ids = chain_entity_id[chain_ids]
    mol.set_annotation('entity_id', entity_ids)
    return entity_ids

def open_structure_rcsb(pdb_code, cache_dir = None):
    import biotite.structure.io.mmtf as mmtf
    import biotite.database.rcsb as rcsb
    
    
    file = mmtf.MMTFFile.read(rcsb.fetch(pdb_code, "mmtf", target_path = cache_dir))
    
    # returns a numpy array stack, where each array in the stack is a model in the 
    # the file. The stack will be of length = 1 if there is only one model in the file
    mol = mmtf.get_structure(file, extra_fields = ["b_factor", "charge"], include_bonds = True) 
    set_atom_entity_id(mol, file)
    return mol, file

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
        
        mol = load(
            pdb_code=pdb_code,
            centre=scene.MN_import_centre, 
            del_solvent=scene.MN_import_del_solvent,
            style=scene.MN_import_style,
            cache_dir=cache_dir, 
            build_assembly = scene.MN_import_build_assembly
        )
        
        bpy.context.view_layer.objects.active = mol
        self.report({'INFO'}, message=f"Imported '{pdb_code}' as {mol.name}")
        
        return {"FINISHED"}

def panel(layout, scene):
    
    layout.label(text = "Download from PDB", icon="IMPORT")
    layout.separator()
    row_import = layout.row()
    row_import.prop(scene, 'MN_pdb_code')
    row_import.operator('mn.import_protein_rcsb', text='Download')
    layout.separator()
    layout.label(text = "Options", icon = "MODIFIER")
    options = layout.column(align = True)
    options.prop(scene, "MN_import_style")
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