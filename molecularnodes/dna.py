import numpy as np
import bpy
from . import obj
from . import coll

bpy.types.Scene.MN_import_oxdna_top = bpy.props.StringProperty(
    name = 'Toplogy File', 
    description = 'File path for the star file to import.', 
    options = {'TEXTEDIT_UPDATE'}, 
    default = '/Users/brady/git/origami/oxdna-viewer/examples/icosahedron/7_icosahedron.json.top', 
    subtype = 'FILE_PATH', 
    maxlen = 0
    )
bpy.types.Scene.MN_import_oxdna_oxdna = bpy.props.StringProperty(
    name = 'oxDNA File', 
    description = 'File path for the oxDNA file to import.', 
    options = {'TEXTEDIT_UPDATE'}, 
    default = '/Users/brady/git/origami/oxdna-viewer/examples/icosahedron/7_icosahedron.json_post_dynamics.oxdna', 
    subtype = 'FILE_PATH', 
    maxlen = 0
    )
bpy.types.Scene.MN_import_oxdna_name = bpy.props.StringProperty(
    name = 'DNA Name', 
    description = 'Name of the created object.', 
    options = {'TEXTEDIT_UPDATE'}, 
    default = 'oxDNA', 
    subtype = 'NONE', 
    maxlen = 0
    )


def read_oxdna(filepath):
    # Open the file and read its contents
    with open(filepath, 'r') as file:
        contents = file.read()

    # Split the contents into lines
    lines = np.array(contents.split('\n'))
    is_meta = np.char.find(lines, '=') > 0

    group_id = np.cumsum(np.append([True], np.diff(is_meta)))
    groups = np.unique(group_id)
    
    frames = []

    for group in groups:
        mask = group == group_id
        if "=" in lines[mask][0]:
            continue
        
        arr = np.loadtxt(lines[mask])
        frames.append(arr)
    
    return np.stack(frames)


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

def add_attributes_to_dna_mol(mol, frame, dna_scale = 0.1):
    attributes = ('base_vector', 'base_normal', 'velocity', 'angular_velocity')
    for i, att in enumerate(attributes):
        col_idx = np.array([3, 4, 5]) + i * 3
        
        try:
            data = frame[:, col_idx]
        except IndexError as e:
            print(f"Unable to get {att} attribute from coordinates. Error: {e}")
            continue
        
        if att != "angular_velocity":
            data *= dna_scale
        
        obj.add_attribute(mol, att, data, type="FLOAT_VECTOR")
    

def load(top, traj, name = 'oxDNA', world_scale = 0.01):
    
    # the scale of the oxDNA files seems to be based on nanometres rather than angstrongs 
    # like most structural biology files, so currently adjusting the world_scale to 
    # compensate
    dna_scale = world_scale * 10
    
    
    frames = read_oxdna(traj)
    topology = read_top(top)
    n_frames = frames.shape[0]
    
    # topology file will return a numeric numpy array, configured as follows:
    #
    # Strand assignment
    # | Base assignment
    # | |  3' Bonded base to the current base (index based on row)
    # | |  |  5' Bonded base to the current base (index based on row)
    # | |  | |
    # S B 3' 5'
    # S B 3' 5'
    # S B 3' 5'
    
    # to get pairs of indices which represent each distinct bond, which are needed for 
    # edge creation in Blender, take each bonded column and create a 'bond' with itself
    idx = np.array(list(range(topology.shape[0])))
    bond_3 = np.vstack((idx, topology[:, 2])).reshape((len(idx), 2))
    bond_5 = np.vstack((idx, topology[:, 3])).reshape((len(idx), 2))
    bonds = np.vstack((bond_3, bond_5))
    
    # drop where either bond is -1 (not bonded) from the bond indices
    mask = bonds == -1
    mask = np.logical_not(mask.any(axis=1)) 

    # creat toplogy object with positions of the first frame, and the bonds from the 
    # topology object
    mol = obj.create_object(
        name=name,
        collection=coll.mn(),
        locations=frames[0][:, 0:3] * dna_scale,
        bonds=bonds[mask, :]
    )
    
    # adding additional toplogy information from the topology and frames objects
    obj.add_attribute(mol, 'res_name', topology[:, 1], "INT")
    obj.add_attribute(mol, 'chain_id', topology[:, 0], "INT")
    add_attributes_to_dna_mol(mol, frames[0], dna_scale=dna_scale)
    
    # if the 'frames' file only contained one timepoint, return the object without creating
    # any kind of collection for storing multiple frames from a trajectory, and a None 
    # object in place of the frames collection
    if n_frames == 1:
        return mol, None
    
    # create a collection to store all of the frame objects that are part of the trajectory
    # they will contain all of the possible attributes which can be interpolated betewen 
    # frames such as position, base_vector, base_normal, velocity, angular_velocity
    collection = coll.frames(f"{name}_frames", parent=coll.data())
    for i, frame in enumerate(frames):
        fill_n = int(np.ceil(np.log10(n_frames)))
        frame_name = f"{name}_frame_{str(i).zfill(fill_n)}"
        frame_mol = obj.create_object(frame_name, collection, frame[:, 0:3] * dna_scale)
        add_attributes_to_dna_mol(frame_mol, frame, dna_scale)
    
    return mol, collection


def panel(layout_function, scene):
    col_main = layout_function.column(heading = "", align = False)
    col_main.label(text = "Import oxDNA File")
    row_import = col_main.row()
    row_import.prop(
        bpy.context.scene, 'MN_import_oxdna_name', 
        text = 'Name', 
        emboss = True
    )
    col_main.prop(bpy.context.scene, 'MN_import_oxdna_top')
    col_main.prop(bpy.context.scene, 'MN_import_oxdna_oxdna')
    row_import.operator('mn.import_oxdna', text = 'Load', icon = 'FILE_TICK')

class MN_OT_Import_Star_File(bpy.types.Operator):
    bl_idname = "mn.import_oxdna"
    bl_label = "Import oxDNA File"
    bl_description = "Will import the given file and toplogy."
    bl_options = {"REGISTER"}

    @classmethod
    def poll(cls, context):
        return True

    def execute(self, context):
        load(
            top=context.scene.MN_import_oxdna_top,
            traj=context.scene.MN_import_oxdna_oxdna, 
            name=context.scene.MN_import_oxdna_name
        )
        return {"FINISHED"}