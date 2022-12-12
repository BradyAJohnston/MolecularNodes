import bpy
import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
from MolecularNodes import load
from MolecularNodes import tools

world_scale = 0.01

def fix_file_cif(file, new_file_name):
    f = open(file_path)
    lines = f.readlines()
    f.close()
    
    for i in range(len(lines)):
        line = lines[i]
        if line.startswith('ATOM'):
            split = line.split()
            split[18] = split[6]
            split[17] = split[16]
            split[16] = split[8]
            new_line = " ".join(split)
            # print(new_line)
            lines[i] = new_line + "\n"
    f = open(new_file_name, "w")
    for line in lines:
        f.write(line)
    f.close()

def add_each_chain(mol, n_chains = 3):
    mol = mol[0]
    
    coll_mn = tools.coll_mn()
    coll_models = bpy.data.collections.get('CellPackComponents')
    if not coll_models:
        coll_models = bpy.data.collections.new('CellPackComponents')
        coll_mn.children.link(coll_models)
    
    chains_unique = np.unique(mol.chain_id)
    chain_number = np.fromiter(
        map(
            lambda x: np.searchsorted(chains_unique, x), 
            mol.chain_id
        ), 
        dtype=int
    )
    
    if len(chains_unique) < n_chains:
        n_chains = len(chains_unique)
    
    for i in range(n_chains):
        mol_sub = mol[chain_number == i]
        mol_name = str(chains_unique[i])
        
        load.create_molecule(
            mol_array = [mol_sub], 
            mol_name = mol_name, 
            collection = coll_models
            )
        

def create_cell(transforms):
    trans_positions = np.array([trans[1] for trans in transforms.values()])
    rot0 = np.array([rot[0][0:1, :][0] for rot in transforms.values()])
    rot1 = np.array([rot[0][1:2, :][0] for rot in transforms.values()])
    rot2 = np.array([rot[0][2:3, :][0] for rot in transforms.values()])

    cell = load.create_object('CellPack', tools.coll_mn(), trans_positions * world_scale)

    rots = [rot1, rot2, rot3]

    for i in range(3):
        rot = "rot" + str(i)
        att = cell.data.attributes.get(rot)
        if not att:
            cell.data.atttributes.new(rot, 'FLOAT_VECTOR', 'POINT')

        att.data.foreach_set('vector', [value for value in rot[i].reshape((len(rot[i]) * 3))])
    

def get_start_end_trans(ops):
    # go through the text and match the starting and ending numbers
    # and return this as a numpy array with a starting column and ending column
    import re
    re_start = re.compile('\d+(?=\-)')
    re_end = re.compile('(?<=\-)\d+')

    arr = np.empty(len(ops) * 2, dtype = int).reshape([len(ops), 2])
    for i in range(len(ops)):
        starts = np.array(re_start.findall(ops[i])).astype(int)
        ends   = np.array(re_end.findall(ops[i])).astype(int)
        starts = starts[0]
        ends = ends[-1]
        arr[i, [0,1]] = [starts, ends]
    return arr

# file_path = "C:\\Users\\BradyJohnston\\Documents\\GitHub\\MycoplasmaGenitalium\\Models\\cellpack_atom_instances_1189_curated\\cellpack_atom_instances_1189_curated.cif"
# new_file = "C:\\Users\\BradyJohnston\\Documents\\GitHub\\MycoplasmaGenitalium\\Models\\cellpack_atom_instances_1189_curated\\cellpack_atom_instances_1189_curated_fixed.cif"

# fix_file_cif(file_path, new_file)

#file = pdbx.PDBxFile.read(new_file)

# mol, file = load.open_structure_local_pdbx(new_file, False)

# add_each_chain(mol, n_chains = 600)

# get_start_end_trans(ops)[:40, :]
# chains_unique = np.unique(chains)

# for i in range(len(ops)):
#     starts = np.array(re_start.findall(ops[i])).astype(int)
#     ends   = np.array(re_end.findall(ops[i])).astype(int)
#     starts = starts[0]
#     ends = ends[0]
#     arr[i, [0,1]] = [starts, ends]