import bpy
import numpy as np
from . import obj, coll


def clear_armature(object):
    for mod in object.modifiers:
        if mod.type == "ARMATURE":
            if mod.object:
                bpy.data.objects.remove(mod.object)
            object.modifiers.remove(mod)


def add_bones(object, name="armature"):
    # creates bones and assigns correct weights

    clear_armature(object)

    bone_positions, bone_weights, chain_ids = get_bone_positions(object)

    armature = create_bones(bone_positions, chain_ids)
    for i in range(bone_weights.shape[1]):
        group = object.vertex_groups.new(name=f"mn_armature_{i}")
        vertex_indices = np.where(bone_weights[:, i] == 1)[0]
        group.add(vertex_indices.tolist(), 1, "ADD")

    object.select_set(True)
    armature.select_set(True)
    bpy.context.view_layer.objects.active = armature
    bpy.ops.object.parent_set(type="ARMATURE")

    bpy.context.view_layer.objects.active = object
    bpy.ops.object.modifier_move_to_index("EXEC_DEFAULT", modifier="Armature", index=0)

    return armature


def get_bone_positions(object):
    positions, atom_name, chain_id, res_id, sec_struct = [
        obj.get_attribute(object, att)
        for att in ["position", "atom_name", "chain_id", "res_id", "sec_struct"]
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
    print("hello world")


def create_bones(positions, chain_ids, name="armature"):
    bpy.ops.object.add(type="ARMATURE", enter_editmode=True)
    object = bpy.context.active_object
    object.name = name
    coll.armature().objects.link(object)
    armature = object.data
    armature.name = f"{name}_frame"
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
        bpy.ops.armature.parent_set(type="CONNECTED")
        for bone in [bone_a, bone_b]:
            armature.edit_bones[bone].select = False
    bpy.ops.object.editmode_toggle()

    return object


class MN_MT_Add_Armature(bpy.types.Operator):
    bl_idname = "mn.add_armature"
    bl_label = "Add Armature"
    bl_description = (
        "Automatically add armature for each amino acid of the structure   "
    )

    def execute(self, context):
        object = context.active_object
        add_bones(bpy.data.objects[object.name], name=object.name)

        return {"FINISHED"}


import bpy
import bmesh
from mathutils import Vector


class CreateVertexGroupAndBone(bpy.types.Operator):
    """Create Vertex Group from Selected Vertices and Corresponding Bone"""

    bl_idname = "object.create_vertex_group_and_bone"
    bl_label = "Create Vertex Group and Bone"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        obj = context.object

        # Ensure we are in edit mode and the object is a mesh
        if obj.mode != "EDIT" or obj.type != "MESH":
            self.report({"ERROR"}, "Must be in Edit Mode on a Mesh Object")
            return {"CANCELLED"}

        # Switch to Object Mode to update selection
        bpy.ops.object.mode_set(mode="OBJECT")

        # Create a new vertex group
        vg = obj.vertex_groups.new(name="My Vertex Group")
        selected_indices = [v.index for v in obj.data.vertices if v.select]
        vg.add(selected_indices, 1.0, "ADD")

        # Calculate median point and principal axis
        bm = bmesh.new()
        bm.from_mesh(obj.data)
        selected_verts = [vert.co for vert in bm.verts if vert.select]
        if not selected_verts:
            self.report({"ERROR"}, "No vertices selected")
            return {"CANCELLED"}
        median_point = sum(selected_verts, Vector()) / len(selected_verts)
        # Simple approximation for bone alignment: Use bounding box
        min_point = Vector(map(min, zip(*selected_verts)))
        max_point = Vector(map(max, zip(*selected_verts)))
        bone_direction = (max_point - min_point).normalized()

        # Check/Create Armature
        bpy.ops.object.mode_set(
            mode="OBJECT"
        )  # Switch to Object Mode to create armature
        armature = None
        for ob in context.scene.objects:
            if ob.type == "ARMATURE" and ob.data.get("My Armature"):
                armature = ob
                break
        if not armature:
            bpy.ops.object.armature_add()
            armature = context.object
            armature.name = "My Armature"
            armature.data.name = "My Armature"

        # Add bone to the armature
        bpy.context.view_layer.objects.active = armature
        bpy.ops.object.mode_set(mode="EDIT")
        bone = armature.data.edit_bones.new(vg.name)
        bone.head = median_point
        bone.tail = median_point + bone_direction  # Simple alignment
        bpy.ops.object.mode_set(mode="OBJECT")

        # Add Armature modifier to mesh and set vertex group
        if "Armature" not in [mod.type for mod in obj.modifiers]:
            mod = obj.modifiers.new("Armature", "ARMATURE")
            mod.object = armature
            mod.use_vertex_groups = True

        # Switch back to Object Mode
        bpy.context.view_layer.objects.active = obj
        bpy.ops.object.mode_set(mode="OBJECT")

        self.report({"INFO"}, "Vertex group and bone created and linked")
        return {"FINISHED"}


# def register():
#     bpy.utils.register_class(CreateVertexGroupAndBone)

# def unregister():
#     bpy.utils.unregister_class(CreateVertexGroupAndBone)

# if __name__ == "__main__":
#     register()

import numpy as np

# Sample data: array of 3D points
points = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]])

# Step 1: Normalize the data
points_mean = points.mean(axis=0)
normalized_points = points - points_mean

# Step 2: Compute the covariance matrix
cov_matrix = np.cov(normalized_points, rowvar=False)

# Step 3: Find the eigenvalues and eigenvectors
eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)

# Step 4: The eigenvector with the highest eigenvalue is the direction
principal_component = eigenvectors[:, np.argmax(eigenvalues)]

print("Principal Component (Overall Trend Direction):", principal_component)
