import bpy
import numpy as np
import MDAnalysis as mda
from . import dict
from .tools import mn_collection
from .open import create_model, add_attribute


def load_trajectory(file_top, 
                    file_traj,
                    world_scale = 0.01, 
                    include_bonds = False, 
                    name = "default"
                    ):
    
    md_start = bpy.context.scene.mol_import_md_frame_start
    md_step =  bpy.context.scene.mol_import_md_frame_step
    md_end =   bpy.context.scene.mol_import_md_frame_end
    
    coll_mn = mn_collection()
    
    # initially load in the trajectory
    univ = mda.Universe(file_top, file_traj)
    
    # TODO allow including of the bonds
    bonds = []
    if include_bonds:
        bonds = []
    
    mol_object = create_model(
        name = name,
        collection = coll_mn, 
        locations = univ.atoms.positions * world_scale, 
        bonds = bonds
    )
    
    is_alpha_carbon = np.isin(
        univ.atoms.ix, 
        univ.select_atoms("name CA").ix
    ).astype(bool)
    
    # add the attributes for the model
    add_attribute(mol_object, 'res_id', univ.atoms.resnums, "INT")
    add_attribute(mol_object, 'b_factor', univ.atoms.tempfactors, "FLOAT")
    add_attribute(mol_object, 'is_alpha_carbon', is_alpha_carbon, "BOOLEAN")
    
    
    
    
    # create the frames of the trajectory in their own collection to be disabled
    coll_frames = bpy.data.collections.new(name + "_frames")
    coll_mn.children.link(coll_frames)
    
    counter = 1
    for ts in univ.trajectory:
        if counter % md_step == 0 and counter >= md_start and counter <= md_end:
            create_model(
                name = name + "_frame_" + str(counter),
                collection = coll_frames, 
                locations = univ.atoms.positions * world_scale
            )
            
            # TODO potential for adding frame-by frame attributes
            # such as energy etc if that information is available
            # will just require a add_attribute() call in this loop
        counter += 1
    
    # disable the frames collection from the viewer
    bpy.context.view_layer.layer_collection.children[coll_mn.name].children[coll_frames.name].exclude = True
    
    return mol_object, coll_frames
    
