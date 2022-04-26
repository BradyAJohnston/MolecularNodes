import bpy
import numpy as np
import re
try:
    import atomium
except:
    print("Atomium Not Installed")

atom_name_dict = {
'C'   : 1,
"C1'" : 2,
'C2'  : 3,
"C2'" : 4,
"C3'" : 5,
'C4'  : 6,
"C4'" : 7,
'C5'  : 8,
"C5'" : 9,
'C6'  : 10,
'C7'  : 11,
'C8'  : 12,
'CA'  : 13,
'CB'  : 14,
'CD'  : 15,
'CD1' : 16,
'CD2' : 17,
'CE'  : 18,
'CE1' : 19,
'CE2' : 20,
'CG'  : 21,
'CG1' : 22,
'CG2' : 23,
'CZ'  : 24,
'N'   : 25,
'N1'  : 26,
'N2'  : 27,
'N3'  : 28,
'N4'  : 29,
'NZ'  : 30,
'O'   : 31,
'O2'  : 32,
"O3'" : 33,
'O4'  : 34,
"O4'" : 35,
"O5'" : 36,
'O6'  : 37,
'OD1' : 38,
'OD2' : 39,
'OE1' : 40,
'OE2' : 41,
'OG'  : 42,
'OG1' : 43,
'OH'  : 44,
'OP1' : 45,
'OP2' : 46,
'OXT' : 47,
'P'   : 48,
'SD'  : 49,
'SG'  : 50
}

element_dict = {
    "H"  : {"atomic_number" : 1, "radii" : 0.53}, 
    "He" : {"atomic_number" : 2, "radii" : 0.31}, 
    "Li" : {"atomic_number" : 3, "radii" : 1.67}, 
    "Be" : {"atomic_number" : 4, "radii" : 1.12}, 
    "B"  : {"atomic_number" : 5, "radii" : 0.87}, 
    "C"  : {"atomic_number" : 6, "radii" : 0.67}, 
    "N"  : {"atomic_number" : 7, "radii" : 0.56}, 
    "O"  : {"atomic_number" : 8, "radii" : 0.48}, 
    "F"  : {"atomic_number" : 9, "radii" : 0.42}, 
    "Ne" : {"atomic_number" : 10, "radii" : 0.38}, 
    "Na" : {"atomic_number" : 11, "radii" : 1.90}, 
    "Mg" : {"atomic_number" : 12, "radii" : 1.45}, 
    "Al" : {"atomic_number" : 13, "radii" : 1.18}, 
    "Si" : {"atomic_number" : 14, "radii" : 1.11}, 
    "P"  : {"atomic_number" : 15, "radii" : 0.98}, 
    "S"  : {"atomic_number" : 16, "radii" : 0.88}, 
    "Cl" : {"atomic_number" : 17, "radii" : 0.79}, 
    "Ar" : {"atomic_number" : 18, "radii" : 0.71}, 
    "K"  : {"atomic_number" : 19, "radii" : 2.43}, 
    "Ca" : {"atomic_number" : 20, "radii" : 1.9}   
}

radii_dict = {
    "H"  : 1.10, 
    "He" : 1.40, 
    "Li" : 1.82, 
    "Be" : 1.53, 
    "B"  : 1.92, 
    "C"  : 1.70, 
    "N"  : 1.55, 
    "O"  : 1.52, 
    "F"  : 1.47, 
    "Ne" : 1.54, 
    "Na" : 2.27, 
    "Mg" : 1.73, 
    "Al" : 1.84, 
    "Si" : 2.10, 
    "P"  : 1.80, 
    "S"  : 1.80, 
    "Cl" : 1.75, 
    "Ar" : 1.88, 
    "K"  : 2.75, 
    "Ca" : 2.31, 
    "Sc" : 2.11, 
    
    # break in the elements, no longer in direct numerical order
    "Ni" : 1.63, 
    "Cu" : 1.40, 
    "Zn" : 1.39
}

AA_dict = {
    # 20 naturally occurring amino acids
    "ALA"  : {"aa_number" : 1, "aa_type" : "polar", "aa_type_no" : 1}, 
    "ARG" : {"aa_number" : 2, "aa_type" : "polar", "aa_type_no" : 1}, 
    "ASN" : {"aa_number" : 3, "aa_type" : "polar", "aa_type_no" : 1}, 
    "ASP" : {"aa_number" : 4, "aa_type" : "polar", "aa_type_no" : 1}, 
    "CYS"  : {"aa_number" : 5, "aa_type" : "polar", "aa_type_no" : 1}, 
    "GLU"  : {"aa_number" : 6, "aa_type" : "polar", "aa_type_no" : 1}, 
    "GLN"  : {"aa_number" : 7, "aa_type" : "polar", "aa_type_no" : 1}, 
    "GLY"  : {"aa_number" : 8, "aa_type" : "polar", "aa_type_no" : 1}, 
    "HIS"  : {"aa_number" : 9, "aa_type" : "polar", "aa_type_no" : 1}, 
    "ILE"  : {"aa_number" : 10, "aa_type" : "polar", "aa_type_no" : 1}, 
    "LEU" : {"aa_number" : 11, "aa_type" : "polar", "aa_type_no" : 1}, 
    "LYS" : {"aa_number" : 12, "aa_type" : "polar", "aa_type_no" : 1}, 
    "MET" : {"aa_number" : 13, "aa_type" : "polar", "aa_type_no" : 1}, 
    "PHE" : {"aa_number" : 14, "aa_type" : "polar", "aa_type_no" : 1}, 
    "PRO" : {"aa_number" : 15, "aa_type" : "polar", "aa_type_no" : 1}, 
    "SER"  : {"aa_number" : 16, "aa_type" : "polar", "aa_type_no" : 1}, 
    "THR"  : {"aa_number" : 17, "aa_type" : "polar", "aa_type_no" : 1}, 
    "TRP" : {"aa_number" : 18, "aa_type" : "polar", "aa_type_no" : 1}, 
    "TYR" : {"aa_number" : 19, "aa_type" : "polar", "aa_type_no" : 1}, 
    "VAL"  : {"aa_number" : 20, "aa_type" : "polar", "aa_type_no" : 1}, 

    # unknown? Came up in one of the structures, haven't looked into it yet
    # TODO look into it!
    "UNK" : {"aa_number" : 21, "aa_type" : "unkown", "aa_type_no" : 1}, 

    ## nucleic acids
    ### DNA
    "DC" : {"aa_number" : 31, "aa_type" : "unkown", "aa_type_no" : 1}, 
    "DG" : {"aa_number" : 32, "aa_type" : "unkown", "aa_type_no" : 1}, 
    "DA" : {"aa_number" : 33, "aa_type" : "unkown", "aa_type_no" : 1}, 
    "DT" : {"aa_number" : 34, "aa_type" : "unkown", "aa_type_no" : 1},

    ### RNA
    "C" : {"aa_number" : 41, "aa_type" : "unkown", "aa_type_no" : 1}, 
    "G" : {"aa_number" : 42, "aa_type" : "unkown", "aa_type_no" : 1}, 
    "A" : {"aa_number" : 43, "aa_type" : "unkown", "aa_type_no" : 1}, 
    "U" : {"aa_number" : 44, "aa_type" : "unkown", "aa_type_no" : 1} 

}

###
# The following variables are passed in when the operators are called, and it is set up properly by Serpens when built
# They are thusly not defined in the actual code, but are included here for reference
###
# create_bonds = True
# pdb_code = "4ozs"
# namometre_scale = 1
# pdb_path = "C:\\Users\\bradyjohnston\\Desktop\\4ozs.pdb"
# connect_cutoff = 0.35

pdb_id = pdb_code
one_nanometre_size_in_metres = nanometre_scale * 0.1

# download the required model
if (fetch_pdb):
    pdb = atomium.fetch(pdb_id)
else: 
    pdb_id = molecule_name
    pdb = atomium.open(pdb_path)

#pdb = atomium.open("C:\\Users\\BradyJohnston\\Desktop\\atp-frames.pdb")

# check the number of models in the file
n_models = len(pdb.models)

first_model = pdb.models[0]
n_atoms = len(first_model.atoms())
all_chains = first_model.chains()

# contains the atom number in the PDB file, acts as the index for everything.
atom_id = []

# contains the XYZ coordinates for the atom
atom_location = []

# contains the character sumbol of the element for the atom, i.e. 'H' for Hydrogen 
# and 'C' for Carbon
atom_element_char = []

# contains the atomic number of the number, has to be matched later against a 
# dictionary as the pip release of atomium doesn't currently have quick-access
# to the atomic number, only the symbol
atom_element_num = []

# contains the name of the atom, which varies depending on where the atom appears
# in the residue, C, C1' etc, with the second list being the numeric encoding of 
# this for use in geometry nodes
atom_name_char = []
atom_name_num = []

# contains the letter code that is the representative of the chain that the 
# residue is part of, with the second list being the numeric encoding of this 
# for use in geometry nodes
atom_chain_char = []
atom_chain_num = []

# contains the 3-letter codes that ID The AA side chain the atom is a part of, 
# or the single or two-letter codes for nucleic acids, with the second list
# being the numeric encoding of this for use in geometry nodes
atom_aa_id_char = []
atom_aa_id_number = []

# contains the integer representing the residue number in the sequence of the 
# protein that the atom is part of
atom_aa_sequence_number = []

# contains the b-factor for each atom, which is a representation of how
# static the atom is when part of the overall structure
atom_b_factor = []

# contains TRUE / FALSE for whether the atom is part of the backbone
atom_is_backbone = []

# contains True / False for whether the atom is part of a sidechain
atom_is_sidechain = []

def try_append(list, value, value_on_fail = 0):
    """
    Tries to append the value to the list, and adds instead the value on fail
    instead if a lookup into one of the dictionaries has failed.
    """
    try:
        list.append(value)
    except:
        list.appen(value_on_fail)


for chain in first_model.chains():
    current_chain = chain.id
    for res in chain.residues():
        current_aa_id_char = res.name
        current_aa_sequence_number = int(re.findall(r"\d+", res.id.split(".")[1])[0])

        for atom in res.atoms():
            # TODO remove the old list appends if everything is working properly
            try_append(atom_id, atom.id)
            # atom_id.append(atom.id)
            try_append(atom_location, atom.location)
            # atom_location.append(atom.location)
            try_append(atom_element_char, atom.element)
            # atom_element_char.append(atom.element)
            try_append(atom_element_num, element_dict[atom.element]["atomic_number"])
            # atom_element_num.append(element_dict[atom.element]["atomic_number"])
            try_append(atom_name_char, atom.name)
            # atom_name_char.append(atom.name)
            try_append(atom_chain_char, current_chain)
            # atom_chain_char.append(current_chain)
            try_append(atom_aa_sequence_number, current_aa_sequence_number)
            # atom_aa_sequence_number.append(current_aa_sequence_number)
            try_append(atom_aa_id_char, current_aa_id_char)
            # atom_aa_id_char.append(current_aa_id_char)
            try_append(atom_b_factor, atom.bvalue)
            # atom_b_factor.append(atom.bvalue)
            try_append(atom_aa_id_number, AA_dict[current_aa_id_char]["aa_number"])
            # atom_aa_id_number.append(AA_dict[current_aa_id_char]["aa_number"])
            try_append(atom_is_backbone, int(atom.is_backbone))
            try_append(atom_is_sidechain, int(atom.is_side_chain))


# turn all of the lists into numpy arrays, so that they can be rearranged based 
# on the atom indices created with np.argsort()

atom_id = np.array(atom_id)
atom_location = np.array(atom_location)
atom_element_char = np.array(atom_element_char)
atom_element_num = np.array(atom_element_num)
atom_name_char = np.array(atom_name_char)
# atom_name_num = np.array(atom_name_num)
atom_chain_char = np.array(atom_chain_char)
# atom_chain_num = np.array(atom_chain_num)
atom_aa_id_char = np.array(atom_aa_id_char)
atom_aa_id_number = np.array(atom_aa_id_number)
atom_aa_sequence_number = np.array(atom_aa_sequence_number)
atom_b_factor = np.array(atom_b_factor)
atom_is_backbone = np.array(atom_is_backbone)
atom_is_sidechain = np.array(atom_is_sidechain)

inds = atom_id.argsort()

# rearrange all of the arrays based on the indices, so that all values 
# match properly, and the atoms are in ascending order
atom_id = atom_id[inds]
atom_location = atom_location[inds]
atom_element_char = atom_element_char[inds]
atom_element_num = atom_element_num[inds]
atom_name_char = atom_name_char[inds]
# atom_name_num = atom_name_num[inds]
atom_chain_char = atom_chain_char[inds]
# atom_chain_num = atom_chain_num[inds]
atom_aa_id_char = atom_aa_id_char[inds]
atom_aa_id_number = atom_aa_id_number[inds]
atom_aa_sequence_number = atom_aa_sequence_number[inds]
atom_b_factor = atom_b_factor[inds]
atom_is_backbone = atom_is_backbone[inds]
atom_is_sidechain = atom_is_sidechain[inds]


unique_chains = np.array(list(set(atom_chain_char)))
chain_inds = unique_chains.argsort()
unique_chains = unique_chains[chain_inds]

atom_chain_num = list(map(lambda x: int(np.where(x == unique_chains)[0]), atom_chain_char))
atom_chain_num = np.array(atom_chain_num)


unique_atoms = np.array(list(set(atom_name_char)))
unique_atoms_inds = unique_atoms.argsort()
unique_atoms = unique_atoms[unique_atoms_inds]


def try_lookup(dict, key, value_on_fail = 0):
    """
    Try looking up the value from the key in the dictionary, 
    otherwise return the value on fail so that things can keep moving
    """
    try:
        return dict[key]
    except:
        return value_on_fail

atom_name_num = list(map(lambda x: int(try_lookup(atom_name_dict, x)), atom_name_char))
atom_name_num = np.array(atom_name_num)


def create_model(name, collection, locations, bonds = [], faces = []):
    """
    Creates a mesh with the given name in the given collection, from the supplied
    values for the locationso of vertices, and if supplied, bonds and faces.
    """
    # create a new mesh
    atom_mesh = bpy.data.meshes.new(name)
    atom_mesh.from_pydata(locations, bonds, faces)
    new_object = bpy.data.objects.new(name, atom_mesh)
    collection.objects.link(new_object)


def create_properties_model(name, collection, prop_x, prop_y, prop_z):
    """
    Creates a mesh that will act as a look up table for properties about the atoms
    in the actual mesh that will be created.
    """
    def get_value(vec, x):
        try:
            return vec[x]
        except:
            return 0

    create_model(
        name = name, 
        collection = collection, 
        locations = list(map(lambda x: [get_value(prop_x, x), get_value(prop_y, x), get_value(prop_z, x)], range(len(atom_aa_sequence_number) - 1))))


def get_bond_list(model, connect_cutoff = 0.35, search_distance = 3):
    """
    For all atoms in the model, search for the nearby atoms given the current 
    distance cutoff, then calculate whether or not they will be bonded to their 
    nearby atoms.

    Returns a list of lists, each with two integers in them, specifying the 
    atoms that are to be bonded.
    """

    mod = model
    mod.optimise_distances()

    for atom in mod.atoms():
        primary_atom = atom
        primary_radius = radii_dict[atom.element]
        nearby_atoms = atom.nearby_atoms(search_distance)
        if atom.element == "H":
            connect_adjust = -0.2
        else:
            connect_adjust = 0

        for secondary_atom in nearby_atoms:
            secondary_radius = radii_dict[secondary_atom.element]
            distance = atom.distance_to(secondary_atom)
            if distance <= ((connect_cutoff + connect_adjust) + (primary_radius + secondary_radius) / 2):
                primary_atom.bond(secondary_atom)


    for atom in mod.atoms():
        if len(atom.bonded_atoms) > 0:
            print(atom.bonded_atoms)

    all_atoms = mod.atoms()
    all_ids = np.array(list(map(lambda x: x.id, all_atoms)))
    inds = all_ids.argsort()
    all_ids = all_ids[inds]

    bond_list = []

    for atom in all_atoms:
        for atom2 in atom.bonded_atoms:
            bond_list.append([int(np.where(atom.id == all_ids)[0]), int(np.where(atom2.id == all_ids)[0])])

    return bond_list


def get_frame_positions(frame):
    """
    Returns a numpy array of all of the atom locations from the given frame. 
    Importantly it orders them according to their atom numbering to sync the frames.
    """
    all_atoms = frame.atoms()
    atom_id = list(map(lambda x: x.id, all_atoms))
    atom_location = list(map(lambda x: x.location, all_atoms))

    atom_id = np.array(atom_id)
    inds = atom_id.argsort()
    atom_id = atom_id[inds]
    atom_location = np.array(atom_location)
    atom_location = atom_location[inds]

    return atom_location



# See if there is a collection called "Molecular Nodes", if so, set it to be the parent
# collection, otherwise create one and link it to the scene collection.

try:
    parent_coll = bpy.data.collections['MolecularNodes']
    parent_coll.name == "MolecularNodes"
    # make the collection active, for creating and disabling new
    bpy.context.view_layer.active_layer_collection = bpy.context.view_layer.layer_collection.children['MolecularNodes']
except:
    parent_coll = bpy.data.collections.new('MolecularNodes')
    bpy.context.scene.collection.children.link(parent_coll)
    # make the collection active, for creating and disabling new
    bpy.context.view_layer.active_layer_collection = bpy.context.view_layer.layer_collection.children['MolecularNodes']



# create new collection that will house the data, link it to the parent collection
col = bpy.data.collections.new(pdb_id)
parent_coll.children.link(col)

col_properties = bpy.data.collections.new(pdb_id + "_properties")
col.children.link(col_properties)

# If create_bonds selected, generate a list of vertex pairs that will be the bonds for the atomic mesh, 
# else return an empty list that will make no edges when passed to create_model()
if create_bonds:
    bonds = get_bond_list(pdb.models[0], connect_cutoff = connect_cutoff)
else:
    bonds = []

# create the first model, that will be the actual atomic model the user will interact with and display
create_model(
    name = pdb_id, 
    collection = col, 
    locations = get_frame_positions(pdb.models[0]) * one_nanometre_size_in_metres, 
    bonds = bonds
    )

# Creat the different models that will encode the various properties into
# the XYZ locations of ther vertices.
create_properties_model(
    name = pdb_id + "_properties_1", 
    collection = col_properties, 
    prop_x = atom_element_num, 
    prop_y = atom_chain_num + 1, # to have the first chain be indexed from 1
    prop_z = atom_name_num
    )
create_properties_model(
    name = pdb_id + "_properties_2", 
    collection = col_properties, 
    prop_x = atom_aa_sequence_number, 
    prop_y = atom_id, 
    prop_z = atom_aa_id_number
    )

create_properties_model(
    name = pdb_id + "_properties_3", 
    collection = col_properties, 
    prop_x = atom_b_factor, 
    prop_y = atom_is_backbone, 
    prop_z = atom_is_sidechain
)

# hide the created properties collection
bpy.context.layer_collection.children[pdb_id].children[pdb_id + '_properties'].exclude = True


if (n_models > 1):
    frames_collection = bpy.data.collections.new(pdb_id + "_frames")
    col.children.link(frames_collection)
    # for each model in the pdb, create a new object and add it to the frames collection
    for frame in pdb.models:
        atom_location = get_frame_positions(frame)
        create_model(
            name = "frame_" + pdb_id, 
            collection = frames_collection, 
            locations = atom_location * one_nanometre_size_in_metres
            )
    
    # hide the created frames collection
    bpy.context.layer_collection.children[pdb_id].children[pdb_id + '_frames'].exclude = True

