import numpy as np
import bpy
import biotite.structure.io.pdbx as pdbx
from . import assembly
from . import obj
from . import load
from . import coll
from . import nodes

bpy.types.Scene.mol_import_cell_pack_path = bpy.props.StringProperty(
    name = 'cellpack_path', 
    description = 'File path for the CellPack file to import.', 
    options = {'TEXTEDIT_UPDATE'}, 
    default = '', 
    subtype = 'FILE_PATH', 
    maxlen = 0
    )
bpy.types.Scene.mol_import_cell_pack_name = bpy.props.StringProperty(
    name = 'cellpack_name', 
    description = 'Name of the created object.', 
    options = {'TEXTEDIT_UPDATE'}, 
    default = 'NewCellPackModel', 
    subtype = 'NONE', 
    maxlen = 0
    )

def load_cellpack(
    file_path, 
    name = 'NewCellPackModel', 
    node_tree = True, 
    world_scale = 0.01
):
    obj_data, coll_cellpack = open_file(file_path, name)
    
    create_cellpack_model(obj_data, coll_cellpack, name = name)

def open_file(file, name = "CellPackModel"):
    file_open = pdbx.PDBxFile.read(file)
    mol = pdbx.get_structure(file_open,  model = 1)
    chain_names = np.unique(mol.chain_id)
    
    # get the transforms and create a data object
    transforms = assembly.cif.CIFAssemblyParser(file_open).get_assemblies()
    obj_data = assembly.mesh.create_data_object(transforms)
    
    coll_cellpack = coll.data("_cellpack")
    
    for chain in chain_names:
        mol_object, coll_frames = load.create_molecule(
            mol_array=mol[mol.chain_id == chain], 
            mol_name=chain, 
            collection=coll_cellpack
            )
        nodes.create_starting_node_tree(mol_object, name = "MOL_cellpack_struc")
    
    return obj_data, coll_cellpack

def create_cellpack_model(obj_data, coll_cellpack, name = "CellPackModel"):
    # create an object with a single vert. This will just the object for instance of the 
    # cellpack data objects
    obj_cellpack = obj.create_object(name = name, collection = coll.mn(), locations=[(0, 0, 0)])
    
    # ensure there is a geometry nodes modifier called 'MolecularNodes' that is created and applied to the object
    node_mod = obj_cellpack.modifiers.get('MolecularNodes')
    if not node_mod:
        node_mod = obj_cellpack.modifiers.new("MolecularNodes", "NODES")

    
    obj_cellpack.modifiers.active = node_mod
    group = nodes.gn_new_group_empty(name = f"MOL_{name}")
    node_mod.node_group = group
    
    # node_mod.node_group = nodes.mol_append_node('MOL_pack_molecules')
    
    node_pack = nodes.add_custom_node_group_to_node(group, 'MOL_pack_molecules')
    node_pack.inputs['Molecules'].default_value = coll_cellpack
    node_pack.inputs['data_object'].default_value = obj_data
    
    link = group.links.new
    link(
        node_pack.outputs[0], 
        group.nodes['Group Output'].inputs[0]
    )


def panel(layout_function, scene):
    col_main = layout_function.column(heading = "", align = False)
    col_main.label(text = "Import CellPack Model")
    row_import = col_main.row()
    row_import.prop(
        bpy.context.scene, 'mol_import_cell_pack_name', 
        text = 'Name', 
        emboss = True
    )
    col_main.prop(
        bpy.context.scene, 'mol_import_cell_pack_path', 
        text = 'CellPack Path (.cif)', 
        emboss = True
    )
    row_import.operator('mol.import_cell_pack', text = 'Load', icon = 'FILE_TICK')

class MOL_OT_Import_Cell_Pack(bpy.types.Operator):
    bl_idname = "mol.import_cell_pack"
    bl_label = "Import CellPack File"
    bl_description = ""
    bl_options = {"REGISTER"}

    @classmethod
    def poll(cls, context):
        return True

    def execute(self, context):
        s = bpy.context.scene
        load_cellpack(
            file_path = s.mol_import_cell_pack_path, 
            name = s.mol_import_cell_pack_name, 
            node_tree = True
        )
        
        return {"FINISHED"}