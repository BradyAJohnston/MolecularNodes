import bpy
import numpy as np
from . import obj

def add_bones(object):
    ## creates bones and assigns correct weights
    
    bone_positions = get_bone_positions(object)
    # weights = get_bone_weights()
    
    bones = create_bones(bone_positions)
    # assign_bone_weights(object, bones, weights)
    
    return bones

def get_bone_positions(object):
    
    positions, atom_name, chain_id, res_id, sec_struct = [
        obj.get_attribute(object, att) 
        for att in ['position', 'atom_name', 'chain_id', 'res_id', 'sec_struct']
        ]
    
    # wherever this is true, a bone will be created and point to the next time it is true
    # we can change this condition. 
    # np.diff(sec_struct) detects changes in secondary structure 
    # atom_name == 2 gets the alpha carbons
    # alpha_carbon_loop = (np.diff(sec_struct) | (sec_struct == 3 & atom_name == 2))
    is_alpha_carbon = atom_name == 2
    idx = np.where(is_alpha_carbon)[0]
    bone_positions = positions[ idx, :]

    print("get_bone_positions")
    return bone_positions

def get_bone_weights(object):
    print("get bone weights")

def create_bones(positions):
    
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

    print("create boens")

def assign_bone_weights(object, bones, weights):
    print("assigning bone weights")
