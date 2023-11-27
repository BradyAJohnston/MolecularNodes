import bpy
import numpy as np
from .blender import obj

def clear_armature(object):
    for mod in object.modifiers:
        if mod.type == "ARMATURE":
            bpy.data.objects.remove(mod.object)
            object.modifiers.remove(mod)

def add_bones(object):
    ## creates bones and assigns correct weights
    
    clear_armature(object)
    
    bone_positions, bone_weights, chain_ids = get_bone_positions(object)
    # weights = get_bone_weights()
    
    armature = create_bones(bone_positions, chain_ids)
    for i in range(bone_weights.shape[1]):
        group = object.vertex_groups.new(name=f'mn_armature_{i}')
        vertex_indices = np.where(bone_weights[:, i] == 1)[0]
        group.add(vertex_indices.tolist(), 1, 'ADD')
    
    object.select_set(True)
    armature.select_set(True)
    bpy.context.view_layer.objects.active = armature
    bpy.ops.object.parent_set(type='ARMATURE')
    
    # for mod in object.modifiers:
    #     if mod.type == "ARMATURE":
            # bpy.ops.object.modifier_move_to_index(modifier=mod.name, index=0)
    bpy.ops.object.modifier_move_to_index('EXEC_DEFAULT', modifier="Armature", index=0)

    
    
    return armature

def get_bone_positions(object):
    
    positions, atom_name, chain_id, res_id, sec_struct = [
        obj.get_attribute(object, att) 
        for att in ['position', 'atom_name', 'chain_id', 'res_id', 'sec_struct']
        ]
    
    is_alpha_carbon = atom_name == 2
    idx = np.where(is_alpha_carbon)[0]
    bone_positions = positions[idx, :]
    bone_positions = np.vstack((bone_positions, positions[-1]))
    group_ids = np.cumsum(is_alpha_carbon)
    groups = np.unique(group_ids)
    bone_weights = np.zeros((len(group_ids), len(groups)))

    for i, unique_id in enumerate(groups):
        bone_weights[:, i] = ((group_ids - 1) == unique_id).astype(int)

    print("get_bone_positions")
    return bone_positions, bone_weights, chain_id[idx]

def get_bone_weights(object):
    print('hello world')

def create_bones(positions, chain_ids):
    
    bpy.ops.object.add(type='ARMATURE', enter_editmode=True)
    object = bpy.context.object
    object.name = 'mn_armature'
    armature = object.data
    armature.name = 'mn_frame'
    arm_name = armature.name
    bones = []
    # add bones
    for i, position in enumerate(positions):
        try:
            pos_a = position
            pos_b = positions[i + 1, :]
        except:
            continue

        bone_name = f"mn_armature_{i}"
        bone = armature.edit_bones.new(bone_name)
        bone.head = pos_a
        bone.tail = pos_b
        bones.append(bone.name)

    armature = bpy.data.armatures[arm_name]
    bones_a = bones.copy()
    bones_b = bones.copy()
    bones_b.pop(0)
    bones = zip(bones_a, bones_b)

    for bone_a, bone_b in bones:
        armature.edit_bones.active = armature.edit_bones[bone_a]
        for bone in [bone_a, bone_b]:
            armature.edit_bones[bone].select = True
        bpy.ops.armature.parent_set(type='CONNECTED')
        for bone in [bone_a, bone_b]:
            armature.edit_bones[bone].select = False

    return object

def assign_bone_weights(object, bones, weights):
    print("assigning bone weights")
