import numpy as np
import bpy
from mathutils import Matrix
from ..blender import  (
    obj, coll
)

def create_data_object(transforms_array, collection = None, name = 'CellPackModel', world_scale = 0.01, fallback=False):
    obj_data = bpy.data.objects.get(name)
    if obj_data and fallback:
        return obj_data
    
    # still requires a unique call TODO: figure out why
    transforms_array = np.unique(transforms_array)
    
    # TODO: this recalculating of chain_ids I don't like, need to figure out a better way
    # to handle this
    chain_ids = np.unique(transforms_array['chain_id'], return_inverse = True)[1] 
    locations = transforms_array['translation'] * world_scale
    
    if not collection:
        collection = coll.data()
    
    obj_data = obj.create_object(locations, collection=collection, name=name)
    obj.add_attribute(obj_data, 'rotation', transforms_array['rotation'], 'QUATERNION', 'POINT')
    obj.add_attribute(obj_data, 'assembly_id', transforms_array['assembly_id'], 'INT', 'POINT')
    obj.add_attribute(obj_data, 'chain_id', chain_ids, 'INT', 'POINT')
    try:
        obj.add_attribute(obj_data, 'transform_id', transforms_array['transform_id'], 'INT', 'POINT')
    except ValueError:
        pass
    
    return obj_data

# data types for the np.array that will store per-chain symmetry operations
dtype = [
    ('assembly_id', int),
    ('transform_id', int),
    ('chain_id',    'U10'),
    ('rotation',  float, 4), # quaternion form
    ('translation', float, 3)
    ]

def array_quaternions_from_dict(transforms_dict):
    n_transforms = 0
    for assembly in transforms_dict.values():
        for transform in assembly:
            n_transforms += len(transform[0])
    
    arr = np.array((n_transforms), dtype=dtype)
    
    transforms = []
    for i, assembly in enumerate(transforms_dict.values()):
        for j, transform in enumerate(assembly):
            chains = transform[0]
            matrix = transform[1]
            arr = np.zeros((len(chains)), dtype = dtype)
            translation, rotation, scale = Matrix(matrix).decompose()
            arr['assembly_id'] = i + 1
            arr['transform_id'] = j
            arr['chain_id'] = chains
            arr['rotation'] = rotation
            arr['translation'] = translation
            transforms.append(arr)
    
    return np.hstack(transforms)
