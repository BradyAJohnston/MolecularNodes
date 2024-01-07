import numpy as np
from pathlib import Path
import bpy

from ..blender import (
    obj, nodes, coll
)
from .load import create_molecule
from .. import utils
from . import parse
from ..color import random_rgb

bpy.types.Scene.mol_import_cell_pack_path = bpy.props.StringProperty(
    name = 'File', 
    description = 'File to import (.cif, .bcif)', 
    subtype = 'FILE_PATH', 
    maxlen = 0
    )
bpy.types.Scene.mol_import_cell_pack_name = bpy.props.StringProperty(
    name = 'Name', 
    description = 'Name of the created object.', 
    default = 'NewCellPackModel', 
    maxlen = 0
    )


def load(
    file_path, 
    name='NewCellPackModel', 
    node_tree=True, 
    world_scale=0.01, 
    fraction: float = 1, 
    instance_nodes=True
    ):
    
    obj_data, chain_collection = read(file_path, name=name, instance_nodes=instance_nodes)
    starting_node_tree(obj_data, chain_collection, name=name, fraction=fraction)
    return obj_data


def read(file, name="NewModel", get_transforms=True, instance_nodes=True):
    if Path(file).suffix in (".bcif", ".bin"):
        data = parse.BCIF(file)
        
        # get transforms and create data / CellPack Object
        if get_transforms:
            obj_data = obj.create_data_object(data.assemblies, name=name, collection=coll.mn())
    else:
        data = parse.PDBX(file, extra_fields = ['label_entity_id'])
        print(f"Loaded model of length {data.n_atoms} with {data.n_chains} chains.")
        
        # get transforms and create data / CellPack Object
        transforms_array = utils.array_quaternions_from_dict(data.assemblies)
        obj_data = obj.create_data_object(transforms_array, name=name, collection=coll.mn())
    
    array = data.structure
    chain_names = np.unique(array.chain_id)
    
    obj_data['chain_id_unique'] = chain_names

    coll_cellpack = coll.cellpack(f"{name}")

    for i, chain in enumerate(chain_names):
        atoms = array[array.chain_id == chain]
        print(chain, atoms.res_name[0])
        
        mol_object, coll_frames = create_molecule(
            array=atoms,
            name=f"{str(i).rjust(4, '0')}_{chain}",
            collection=coll_cellpack
            )

        colors = np.tile(random_rgb(i), (len(atoms), 1))

        obj.add_attribute(mol_object, name="Color", data=colors, type="FLOAT_COLOR", overwrite=True)
        if instance_nodes:
            nodes.create_starting_node_tree(mol_object, name = f"MN_pack_instance_{name}", set_color=False)

    return obj_data, coll_cellpack


def starting_node_tree(ensemble, coll_cellpack, name = "CellPackModel", fraction: float = 1.0, fallback=False):
    # ensure there is a geometry nodes modifier called 'MolecularNodes' that is created and applied to the object
    node_mod = ensemble.modifiers.get('MolecularNodes')
    if not node_mod:
        node_mod = ensemble.modifiers.new("MolecularNodes", "NODES")

    ensemble.modifiers.active = node_mod
    
    group = nodes.new_group(name = f"MN_cellpack_{name}", fallback=False)
    node_mod.node_group = group
    
    node_pack = nodes.add_custom(group, 'MN_pack_instances', location=[-100,0])
    node_pack.inputs['Collection'].default_value = coll_cellpack
    node_pack.inputs['Fraction'].default_value = fraction
    
    link = group.links.new
    link(nodes.get_input(group).outputs[0], node_pack.inputs[0])
    link(node_pack.outputs[0], nodes.get_output(group).inputs[0])


class MN_OT_Import_Cell_Pack(bpy.types.Operator):
    bl_idname = "mol.import_cell_pack"
    bl_label = "Load"
    bl_description = ""
    bl_options = {"REGISTER"}

    @classmethod
    def poll(cls, context):
        return True

    def execute(self, context):
        s = context.scene
        load(
            file_path = s.mol_import_cell_pack_path, 
            name = s.mol_import_cell_pack_name, 
            node_tree = True
        )
        return {"FINISHED"}

def panel(layout, scene):
    layout.label(text = "Load CellPack Model", icon='FILE_TICK')
    layout.separator()
    row_import = layout.row()
    row_import.prop(scene, 'mol_import_cell_pack_name')
    layout.prop(scene, 'mol_import_cell_pack_path')
    row_import.operator('mol.import_cell_pack')