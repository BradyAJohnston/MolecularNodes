import numpy as np
import bpy
import biotite.structure.io.pdbx as pdbx
from . import assembly
from . import obj
from . import load
from . import coll
from . import nodes
from .color import random_rgb
from pathlib import Path
from . import bcif

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


def load_cellpack(file_path, 
                  name = 'NewCellPackModel', 
                  node_tree = True, 
                  world_scale = 0.01
                  ):
    obj_data, coll_cellpack = open_file(file_path, name=name)
    
    starting_node_tree(obj_data, coll_cellpack, name = name)


def open_file(file, name="NewModel", get_transforms=True):
    print("openfile",file)
    if Path(file).suffix in (".bcif", ".bin"):
        mol, transforms = bcif.parse(file)
    else:
        file_open = pdbx.PDBxFile.read(file)
        print("file_open ok")
        mol = pdbx.get_structure(file_open,  model=1, extra_fields=['label_entity_id'])
        print("loaded mol", len(mol))
        transforms = assembly.cif.CIFAssemblyParser(file_open).get_assemblies()
        print("loaded transforms", len(transforms))
    chain_names = np.unique(mol.chain_id)
    # get the transforms and create a data object
    if get_transforms:
        obj_data = assembly.mesh.create_data_object(transforms, name=name)

    coll_cellpack = coll.cellpack(f"{name}")

    for i, chain in enumerate(chain_names):
        atoms = mol[mol.chain_id == chain]
        print(chain, atoms.res_name[0])
        # if atoms.res_name[0] == 'LIP': 
        #     continue
        mol_object, coll_frames = load.create_molecule(
            MN_array=atoms,
            MN_name=f"{str(i).rjust(4, '0')}_{chain}",
            collection=coll_cellpack
            )

        colors = np.tile(random_rgb(), (len(atoms), 1))

        obj.add_attribute(mol_object, name="Color", data=colors, type="FLOAT_COLOR", overwrite=True)
        nodes.create_starting_node_tree(mol_object, name = f"MN_pack_instance_{name}", set_color=False)

    return obj_data, coll_cellpack

def starting_node_tree(obj_data, coll_cellpack, name = "CellPackModel", fraction: float = 1.0, fallback=False):
    # create an object with a single vert. This will just the object for instance of the 
    # cellpack data objects
    
    # ensure there is a geometry nodes modifier called 'MolecularNodes' that is created and applied to the object
    node_mod = obj_data.modifiers.get('MolecularNodes')
    if not node_mod:
        node_mod = obj_data.modifiers.new("MolecularNodes", "NODES")

    obj_data.modifiers.active = node_mod
    
    group = nodes.gn_new_group_empty(name = f"MN_cellpack_{name}", fallback=False)
    node_mod.node_group = group
    
    node_pack = nodes.add_custom_node_group_to_node(group, 'MN_pack_instances')
    node_pack.inputs['Collection'].default_value = coll_cellpack
    node_pack.inputs['Fraction'].default_value = fraction
    
    link = group.links.new
    link(
        group.nodes['Group Input'].outputs[0], 
        node_pack.inputs[0]
    )
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

class MN_OT_Import_Cell_Pack(bpy.types.Operator):
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