import numpy as np
from . import obj
from . import coll

# filepath = '/Users/brady/git/origami/oxdna-viewer/examples/icosahedron/7_icosahedron.json.oxdna'
filepath = '/Users/brady/git/origami/oxdna-viewer/examples/icosahedron/7_icosahedron.json_post_dynamics.oxdna'

def read_oxdna(file):
    # Open the file and read its contents
    with open(filepath, 'r') as file:
        contents = file.read()

    # Split the contents into lines
    lines = np.array(contents.split('\n'))
    is_meta = np.char.find(lines, '=') > 0

    group_id = np.cumsum(np.append(np.array([True]), np.diff(is_meta)))
    groups = np.unique(group_id)
    
    frames = []

    for group in groups:
        mask = group == group_id
        if "=" in lines[mask][0]:
            continue
        
        arr = np.loadtxt(lines[mask])
        frames.append(arr)
    
    return np.stack(frames)

print(read_oxdna(filepath)[0].shape)

top_file_path = '/Users/brady/git/origami/oxdna-viewer/examples/icosahedron/7_icosahedron.json.top'

def read_top(filepath):
    with open(filepath, 'r') as file:
        contents = file.read()

    lines = np.array(contents.split('\n'))
    meta = lines[0]
    arr = np.loadtxt(lines[1:], dtype=str)
    # convert the columns to numeric
    arr_numeric = np.zeros(arr.shape, dtype = int)
    arr_numeric[:, (0, 2, 3)] = arr[:, (0, 2, 3)].astype(int)
    arr_numeric[:, 1] = np.unique(arr[:, 1], return_inverse=True)[1] + 30 # add offset for numeric rep of bases
    # arr_numeric[arr_numeric == -1] = arr_numeric.shape[0]
    
    return arr_numeric

def load(top, traj, world_scale = 0.01):
    
    dna_scale = world_scale * 10
    
    frames = read_oxdna(traj)[0]
    topology = read_top(top)
    idx = np.array(list(range(topology.shape[0])))
    bond_3 = np.vstack((idx, topology[:, 2])).reshape((len(idx), 2))
    bond_5 = np.vstack((idx, topology[:, 3])).reshape((len(idx), 2))
    bonds = np.vstack((bond_3, bond_5))
    mask = bonds == -1
    mask = np.logical_not(mask.any(axis=1))  # Invert the boolean value

    mol = obj.create_object(
        name='oxdna',
        collection=coll.mn(),
        locations=frames[:, 0:3] * dna_scale,
        bonds=bonds[mask, :]
    )
    
    obj.add_attribute(mol, 'res_name', topology[:, 1], "INT")
    obj.add_attribute(mol, 'strand_idx', topology[:, 0], "INT")
    
    attributes = ('base_vector', 'base_normal', 'velocity', 'angular_velocity')
    
    for i, att in enumerate(attributes):
        col_idx = np.array([3, 4, 5]) + i * 3
        
        data = frames[:, col_idx]
        if "velocity" in att:
            data *= dna_scale
        
        obj.add_attribute(mol, att, data, type="FLOAT_VECTOR")

def yay():
    load(top_file_path, filepath)