import numpy as np
import biotite.structure as struc

cols = (
    "group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id",
    "label_comp_id", "label_asym_id", "label_entity_id", "label_seq_id",
    "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy",
    "B_iso_or_equiv", "pdbx_formal_charge", "auth_seq_id", "auth_comp_id",
    "auth_asym_id", "auth_atom_id", "pdbx_PDB_model_num"
)

def rotation_from_matrix(matrix):
    from scipy.spatial.transform import Rotation as R
    
    # calculate the euler rotation from the rotation matrix
    # Blender is 'xyz' euler rotations. Internally they use matrices / quaternions, but
    # current interfaces for geometry nodes are just eulers
    
    rotation = R.from_matrix(matrix).as_euler('xyz')
    
    return rotation

def parse(file):
    dat = read(file)
    arr = arr_from_dat(dat)
    syms = get_syms(dat)
    return arr, syms

def read(file):
    from mmcif.io.BinaryCifReader import BinaryCifReader
    reader = BinaryCifReader()
    dat = reader.deserialize(file)
    return dat

def arr_from_dat(dat):
    atom_site = np.array(dat[0].getObj('atom_site'))
    
    n_atom = len(atom_site)
    arr = struc.AtomArray(n_atom)

    def get_col(col_names):
        return atom_site[:, np.where(np.isin(cols, col_names))[0]].reshape(n_atom)

    coords = np.hstack([
            get_col(f'Cartn_{axis}').astype(float).reshape((n_atom, 1)) for axis in "xyz"
        ])
    
    arr.coord = coords
    arr.set_annotation('chain_id', get_col('label_asym_id').astype('<U4'))
    arr.atom_name = get_col('label_atom_id')
    arr.element = get_col('type_symbol')
    arr.res_name = get_col('label_comp_id')
    
    return arr

def expand_grid(x, y):
    grid_x, grid_y = np.meshgrid(x, y)
    return np.column_stack((grid_x.ravel(), grid_y.ravel()))

def get_syms(dat):
    """Convert symmetry operations to numpy array from data object
    """
    assembly_gen = np.array(dat[0].getObj('pdbx_struct_assembly_gen'))
    ops        = np.array(dat[0].getObj('pdbx_struct_oper_list'))
    
    assemblies = np.char.split(np.char.strip(assembly_gen[:, 1], ('()')), '-')
    
    dtype = [
        ('assembly_id', int),
        ('chain_id',    'U10'),
        ('trans_id', int),
        ('rotation',    float, 3),
        ('translation', float, 3)
    ]
    lists = []
    masks = []
    for i, assembly in enumerate(assemblies):
        start, end = np.array(assembly).astype(int)
        mask = np.array(range(start, end + 1)) - 1
        masks.append(mask)
        
        chain_ids = np.array(assembly_gen[i, 2].split(','))
        lists.append(expand_grid(mask, chain_ids))
    
    all_ops = np.row_stack(lists)
    operations = np.zeros(len(all_ops), dtype = dtype)
    operations['trans_id'] = all_ops[:, 0]
    operations['chain_id'] = all_ops[:, 1]
    
    for i, sym in enumerate(ops):
        mat = sym[1:13].reshape((3, 4), order = 'F')
        rotation = rotation_from_matrix(mat[:3, :3])
        translation = mat[:3, 3]
        mask = operations['trans_id'] == i
        operations['rotation'][mask] = rotation
        operations['translation'][mask] = translation
    
    return operations
