import bpy
from . import obj
from . import load
from . import nodes
import requests
import io

bpy.types.Scene.MN_esmfold_sequence = bpy.props.StringProperty(
    name = 'amino_acid_sequence', 
    description = 'Amino acid sequence of the structure to open', 
    options = {'TEXTEDIT_UPDATE'}, 
    default = '', 
    subtype = 'FILE_PATH', 
    maxlen = 0
    )
bpy.types.Scene.MN_esmfold_name = bpy.props.StringProperty(
    name = 'MN_name', 
    description = 'Name of the molecule on import', 
    options = {'TEXTEDIT_UPDATE'}, 
    default = 'NewMolecule', 
    subtype = 'NONE', 
    maxlen = 0
    )

def molecule_esmfold(
    amino_acid_sequence,               
    MN_name = "Name",                   
    center_molecule = False,               
    del_solvent = True,               
    include_bonds = True,   
    starting_style = 'atoms',               
    setup_nodes = True              
    ):
    mol, file = open_structure_esm_fold(
        amino_acid_sequence = amino_acid_sequence, 
        include_bonds=include_bonds
        )
    
    MN_object, coll_frames = load.create_molecule(
        MN_array = mol,
        MN_name = MN_name,
        file = file,
        calculate_ss = True,
        center_molecule = center_molecule,
        del_solvent = del_solvent, 
        include_bonds = True
        )
    
    if setup_nodes:
        nodes.create_starting_node_tree(
            obj = MN_object, 
            coll_frames=coll_frames, 
            starting_style = starting_style
            )    
    return MN_object

def open_structure_esm_fold(amino_acid_sequence, include_bonds=True):
    import biotite.structure.io.pdb as pdb
    
    
    
    # output_of_subprocess = subprocess.Popen([
    # 'curl',
    # '-X',
    # 'POST',
    # '--data',
    # amino_acid_sequence,
    # 'https://api.esmatlas.com/foldSequence/v1/pdb/'
    # ], stdout=subprocess.PIPE)

    # (esm_folded_pdb_str, err) = output_of_subprocess.communicate()

    r = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', data=amino_acid_sequence)
    
    if r.ok:
    
        esm_folded_pdb_str = r.text
        with io.StringIO() as f:
            f.write(esm_folded_pdb_str)
            f.seek(0)
            file = pdb.PDBFile.read(f)
    
        # returns a numpy array stack, where each array in the stack is a model in the 
        # the file. The stack will be of length = 1 if there is only one model in the file
        mol = pdb.get_structure(file, extra_fields = ['b_factor', 'charge'], include_bonds = include_bonds)
        return mol, file

    else:
        raise ValueError(f'ESMFold returned an error for the amino acid sequence input. This is the error message: {r.text}')

# operator that calls the function to import the structure from ESMFold
class MN_OT_Import_Protein_ESMFold(bpy.types.Operator):
    bl_idname = "mn.import_protein_esmfold"
    bl_label = "import_protein_esmfold"
    bl_description = "Generate structure from ESMFold"
    bl_options = {"REGISTER", "UNDO"}

    @classmethod
    def poll(cls, context):
        return not False

    def execute(self, context):
        amino_acid_sequence = bpy.context.scene.MN_esmfold_sequence
        
        MN_object = molecule_esmfold(
            amino_acid_sequence=amino_acid_sequence, 
            MN_name=bpy.context.scene.MN_esmfold_name,
            include_bonds=bpy.context.scene.MN_import_include_bonds, 
            center_molecule=bpy.context.scene.MN_import_center, 
            del_solvent=bpy.context.scene.MN_import_del_solvent, 
            starting_style=bpy.context.scene.MN_import_default_style, 
            setup_nodes=True
            )
        
        # return the good news!
        bpy.context.view_layer.objects.active = MN_object
        self.report({'INFO'}, message=f"Generated protein '{amino_acid_sequence}' as {MN_object.name}")
        return {"FINISHED"}

    def invoke(self, context, event):
        return self.execute(context)

# UI panel definition

def panel(layout_function, ):
    col_main = layout_function.column(heading = '', align = False)
    col_main.alert = False
    col_main.enabled = True
    col_main.active = True
    col_main.label(text = "Generate Structure from ESMFold")
    row_name = col_main.row(align = False)
    row_name.prop(bpy.context.scene, 'MN_esmfold_name', 
                    text = "Name", icon_value = 0, emboss = True)
    row_name.operator('mn.import_protein_esmfold', text='Generate', icon='IMPORT')
    
    row_seq = col_main.row()
    row_seq.prop(
        bpy.context.scene, 'MN_esmfold_sequence', 
        text = "Sequence", 
        icon_value = 0, 
        emboss = True
    )