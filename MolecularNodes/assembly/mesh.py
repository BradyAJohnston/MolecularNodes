import numpy as np
import bpy
from .. import obj
from .. import coll

def create_data_object(transforms_dict, name = 'DataObject', world_scale = 0.01):
    obj_data = bpy.data.objects.get(name)
    if obj_data:
        return obj_data
    
    transforms_array = get_transforms_from_dict(transforms_dict)
    chain_ids = np.unique(transforms_array['chain_id'], return_inverse = True)[1] 
    locations = transforms_array['translation'] * world_scale
    
    obj_data = obj.create_object(name, coll.data(), locations)
    obj.add_attribute(obj_data, 'assembly_rotation', transforms_array['rotation'], 'FLOAT_VECTOR', 'POINT')
    obj.add_attribute(obj_data, 'assembly_id', transforms_array['assembly_id'], 'INT', 'POINT')
    obj.add_attribute(obj_data, 'chain_id', chain_ids, 'INT', 'POINT')
    
    return obj_data

# data types for the np.array that will store per-chain symmetry operations
dtype = [
    ('assembly_id', int),
    ('chain_id',    'U10'),
    ('rotation',    float, 3),
    ('translation', float, 3)
    ]

def get_transforms_from_dict(transforms_dict):
    
    results = None
    for assembly_id, transforms in transforms_dict.items():
        transforms = transforms_from_assemblies(transforms, index = int(assembly_id))
        if results is None:
            results = transforms
        else:
            results = np.concatenate((results, transforms), axis = 0)
    
    # currently also using unique to remove duplicates, TODO: look into duplicates
    results = np.unique(results)
    
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
        current_transform += n_chains
        transforms = transform_chains(assembly, index = index)
        results[mask] = transforms
    
    return results

def transform_chains(assembly, index = 0):
    
    chains = assembly[0]
    rotation_matrix = np.array(assembly[1]).reshape((3, 3))
    rotation_euler = rotation_from_matrix(rotation_matrix)
    translation_matrix = np.array(assembly[2])
    
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
    import warnings
    
    # calculate the euler rotation from the rotation matrix
    # Blender is 'xyz' euler rotations. Internally they use matrices / quaternions, but
    # current interfaces for geometry nodes are just eulers
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return R.from_matrix(matrix).as_euler('xyz')