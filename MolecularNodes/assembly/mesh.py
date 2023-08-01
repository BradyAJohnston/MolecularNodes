import numpy as np
import bpy

dtype = [('assembly_id', int),
         ('chain_id', 'U10'),
         ('rotation', float, 3),
         ('translation', float, 3)]

def get_transforms_from_dict(transforms_dict):
    
    results = None
    for assembly_id, transforms in transforms_dict.items():
        transforms = transforms_from_assemblies(transforms, index = int(assembly_id))
        if results is None:
            results = transforms
        else:
            results = np.concatenate((results, transforms), axis = 0)
    
    results = results[results['chain_id'] != ''] # filter out where chain is '0', TODO look into why this happens
    results = results.copy(order = 'c') # blender only likes numpy arrays in this order
    
    return results

def transforms_from_assemblies(assembly_list, index = 0):
    n_transforms = 0
    for assembly in assembly_list:
        n_transforms += len(assembly[0])
    
    results = np.zeros((n_transforms), dtype = dtype)
    
    current_transform = 0
    for i, assembly in enumerate(assembly_list):
        n_chains = len(assembly[0])
        
        mask = np.array(range(n_chains))
        mask += current_transform
        
        transforms = transform_chains(assembly, index = index)
        results[mask] = transforms
    
    return results

def transform_chains(assembly, index = 0):
    
    chains = assembly[0]
    rotation_matrix = assembly[1]
    rotation_euler = rotation_from_matrix(rotation_matrix)
    translation_matrix = assembly[2]
    
    n = len(chains)
    result = np.zeros((n), dtype = dtype)
    
    for i, chain in enumerate(chains):
        result[i]['assembly_id'] = index
        result[i]['chain_id'] = chain
        result[i]['rotation'] = rotation_euler
        result[i]['translation'] = translation_matrix
        
    return result

def rotation_from_matrix(matrix):
    from scipy.spatial.transform import Rotation as R
    
    # calculate the euler rotation from the rotation matrix
    # Blender is 'xyz' euler rotations. Internally they use matrices / quaternions, but
    # current interfaces for geometry nodes are just eulers
    
    rotation = R.from_matrix(matrix).as_euler('xyz')
    
    return rotation