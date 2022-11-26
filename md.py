import bpy
import numpy as np
from . import data
from .tools import coll_mn
from .load import create_object, add_attribute
import warnings


def load_trajectory(file_top, 
                    file_traj,
                    world_scale = 0.01, 
                    include_bonds = False, 
                    name = "default"
                    ):
    
    import MDAnalysis as mda
    import MDAnalysis.transformations as trans
    
    md_start = bpy.context.scene.mol_import_md_frame_start
    md_step =  bpy.context.scene.mol_import_md_frame_step
    md_end =   bpy.context.scene.mol_import_md_frame_end
    
    # initially load in the trajectory
    univ = mda.Universe(file_top, file_traj)
    
    bonds = []
    if include_bonds:
        # TODO: support bonds in MD import
        # potentially have to recalculate, MDAnalysis has some support for bond calculation
        bonds = []
    
    
    # fix issues with proteins going outside of the bounding box of the simulation
    # TODO: ask about why this might not be working
    center_on_protein = False
    
    if center_on_protein:
        try:
            protein = univ.select_atoms('protein')
            for ts in univ.trajectory:
                protein.unwrap(compound = 'fragments')
                protein_center = protein.center_of_mass(wrap = True)
                dim = ts.triclinic_dimensions
                box_center = np.sum(dim, axis = 0) / 2
                univ.atoms.translate(box_center - protein_center)
        except:
            warnings.warn("Unable to center the protein during import.")

    # wrap the solvent around as well to better wrap around the protein
    wrap_molecules = False
    if wrap_molecules:
        try:
            not_protein = univ.select_atoms('not protein')
            for ts in univ.trajectory:
                not_protein.wrap(compound = 'residues')
        except:
            warnings.warn("Unable to wrap molecules on import.")



    # create the initial model
    mol_object = create_object(
        name = name,
        collection = coll_mn(), 
        locations = univ.atoms.positions * world_scale, 
        bonds = bonds
    )
    
    
    # add the attributes for the model
    
    
    
    # atomic_number
    try:
        try:
            atomic_name = univ.atoms.elements
        except:
            atomic_name = np.array(list(map(
                lambda x: mda.topology.guessers.guess_atom_element(x), 
                univ.atoms.names
                )))
        atomic_number = np.array(list(map(
            lambda x: data.elements.get(x, {"atomic_number": -1}).get("atomic_number"), 
            np.char.title(atomic_name)
        )))
        add_attribute(mol_object, 'atomic_number', atomic_number, 'INT')
    except:
        warnings.warn("Unable to add atomic numbers")
    
    # vdw radii
    try:
        vdw_radii = np.array(list(map(
            lambda x: mda.topology.tables.vdwradii.get(x, 1), 
            np.char.upper(atomic_name)
        )))
        
        vdw_radii = vdw_radii * world_scale
        
        add_attribute(mol_object, 'vdw_radii', vdw_radii, 'FLOAT')
    except:
        warnings.warn("Unable to add vdw_radii attribute.")
    
    # residue ids (numeric)
    try:
        add_attribute(mol_object, 'res_id', univ.atoms.resnums, "INT")
    except:
        warnings.warn("Unable to add residue IDs")
    
    ### residue names converted to integers in alphabetical order
    try:
        res_names =  np.array(list(map(lambda x: x[0: 3], univ.atoms.resnames)))
        res_numbers = np.array(list(map(lambda x: data.residues.get(x, {'res_name_num': 0}).get('res_name_num'), res_names)))
        add_attribute(mol_object, 'res_name', res_numbers, "INT")
    except:
        warnings.warn("Unable to add residue names")
    
    # TODO: find better way to handle univ that is missing components
    # this works currently but its ugly, need better way to handle missing
    # components and warng the user rather than a massive error
    try:
        add_attribute(mol_object, 'b_factor', univ.atoms.tempfactors, "FLOAT")
    except:
        warnings.warn("Unable to add b_factor (may not be supported in this structure")
    
    # chain_id
    try:
        chain_id = univ.atoms.chainIDs
        chain_id_unique = np.unique(chain_id)
        
        chain_id_num = np.array(list(map(lambda x: np.where(x == chain_id_unique)[0][0], chain_id)))
        add_attribute(mol_object, 'chain_id', chain_id_num, "INT")
    except:
        warnings.warn("Unable to add Chain IDs")
        
    
    # is_backbone
    try:
        is_backbone = np.isin(univ.atoms.ix, univ.select_atoms("backbone or nucleicbackbone").ix).astype(bool)
        add_attribute(mol_object, 'is_backbone', is_backbone, "BOOLEAN")
    except:
        warnings.warn("Unable to add is_backbone attribute.")
    
    try:
        is_alpha_carbon = np.isin(univ.atoms.ix, univ.select_atoms("name CA").ix).astype(bool)
        add_attribute(mol_object, 'is_alpha_carbon', is_alpha_carbon, "BOOLEAN")
    except:
        warnings.warn("Unable to add is_ca attribute.")
    
    
    
    
    # create the frames of the trajectory in their own collection to be disabled
    coll_frames = bpy.data.collections.new(name + "_frames")
    coll_mn().children.link(coll_frames)
    
    for ts in univ.trajectory[md_start:md_end:md_step]:
        create_object(
            name = name + "_frame_" + str(ts.frame),
            collection = coll_frames, 
            locations = univ.atoms.positions * world_scale
        )

        # TODO potential for adding frame-by frame attributes
        # such as energy etc if that information is available
        # will just require a add_attribute() call in this loop
    
    # disable the frames collection from the viewer
    bpy.context.view_layer.layer_collection.children[coll_mn().name].children[coll_frames.name].exclude = True
    
    return mol_object, coll_frames
    
