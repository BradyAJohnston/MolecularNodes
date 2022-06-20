import bpy
import numpy as np
import re
import sys
import os
import site


def verify_user_sitepackages():
    usersitepackagespath = site.getusersitepackages()

    if os.path.exists(usersitepackagespath) and usersitepackagespath not in sys.path:
        sys.path.append(usersitepackagespath)


verify_user_sitepackages()

try:
    import atomium
except:
    print("Atomium Not Installed")

atom_name_dict = {
    'C': 1,
    "C1'": 2,
    'C2': 3,
    "C2'": 4,
    "C3'": 5,
    'C4': 6,
    "C4'": 7,
    'C5': 8,
    "C5'": 9,
    'C6': 10,
    'C7': 11,
    'C8': 12,
    'CA': 13,
    'CB': 14,
    'CD': 15,
    'CD1': 16,
    'CD2': 17,
    'CE': 18,
    'CE1': 19,
    'CE2': 20,
    'CG': 21,
    'CG1': 22,
    'CG2': 23,
    'CZ': 24,
    'N': 25,
    'N1': 26,
    'N2': 27,
    'N3': 28,
    'N4': 29,
    'NZ': 30,
    'O': 31,
    'O2': 32,
    "O3'": 33,
    'O4': 34,
    "O4'": 35,
    "O5'": 36,
    'O6': 37,
    'OD1': 38,
    'OD2': 39,
    'OE1': 40,
    'OE2': 41,
    'OG': 42,
    'OG1': 43,
    'OH': 44,
    'OP1': 45,
    'OP2': 46,
    'OXT': 47,
    'P': 48,
    'SD': 49,
    'SG': 50
}

element_dict = {
    "H": {"atomic_number": 1,  "radii": 1.10},
    "He": {"atomic_number": 2,  "radii": 1.40},
    "Li": {"atomic_number": 3,  "radii": 1.82},
    "Be": {"atomic_number": 4,  "radii": 1.53},
    "B": {"atomic_number": 5,  "radii": 1.92},
    "C": {"atomic_number": 6,  "radii": 1.70},
    "N": {"atomic_number": 7,  "radii": 1.55},
    "O": {"atomic_number": 8,  "radii": 1.52},
    "F": {"atomic_number": 9,  "radii": 1.47},
    "Ne": {"atomic_number": 10, "radii": 1.54},
    "Na": {"atomic_number": 11, "radii": 2.27},
    "Mg": {"atomic_number": 12, "radii": 1.73},
    "Al": {"atomic_number": 13, "radii": 1.84},
    "Si": {"atomic_number": 14, "radii": 2.10},
    "P": {"atomic_number": 15, "radii": 1.80},
    "S": {"atomic_number": 16, "radii": 1.80},
    "Cl": {"atomic_number": 17, "radii": 1.75},
    "Ar": {"atomic_number": 18, "radii": 1.88},
    "K": {"atomic_number": 19, "radii": 2.75},
    "Ca": {"atomic_number": 20, "radii": 2.31}
}

radii_dict = {
    "H": 1.10,
    "He": 1.40,
    "Li": 1.82,
    "Be": 1.53,
    "B": 1.92,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "Ne": 1.54,
    "Na": 2.27,
    "Mg": 1.73,
    "Al": 1.84,
    "Si": 2.10,
    "P": 1.80,
    "S": 1.80,
    "Cl": 1.75,
    "Ar": 1.88,
    "K": 2.75,
    "Ca": 2.31,
    "Sc": 2.11,

    # break in the elements, no longer in direct numerical order
    "Ni": 1.63,
    "Cu": 1.40,
    "Zn": 1.39
}

AA_dict = {
    # 20 naturally occurring amino acids
    "ALA": {"aa_number": 1, "aa_type": "polar", "aa_type_no": 1},
    "ARG": {"aa_number": 2, "aa_type": "polar", "aa_type_no": 1},
    "ASN": {"aa_number": 3, "aa_type": "polar", "aa_type_no": 1},
    "ASP": {"aa_number": 4, "aa_type": "polar", "aa_type_no": 1},
    "CYS": {"aa_number": 5, "aa_type": "polar", "aa_type_no": 1},
    "GLU": {"aa_number": 6, "aa_type": "polar", "aa_type_no": 1},
    "GLN": {"aa_number": 7, "aa_type": "polar", "aa_type_no": 1},
    "GLY": {"aa_number": 8, "aa_type": "polar", "aa_type_no": 1},
    "HIS": {"aa_number": 9, "aa_type": "polar", "aa_type_no": 1},
    "ILE": {"aa_number": 10, "aa_type": "polar", "aa_type_no": 1},
    "LEU": {"aa_number": 11, "aa_type": "polar", "aa_type_no": 1},
    "LYS": {"aa_number": 12, "aa_type": "polar", "aa_type_no": 1},
    "MET": {"aa_number": 13, "aa_type": "polar", "aa_type_no": 1},
    "PHE": {"aa_number": 14, "aa_type": "polar", "aa_type_no": 1},
    "PRO": {"aa_number": 15, "aa_type": "polar", "aa_type_no": 1},
    "SER": {"aa_number": 16, "aa_type": "polar", "aa_type_no": 1},
    "THR": {"aa_number": 17, "aa_type": "polar", "aa_type_no": 1},
    "TRP": {"aa_number": 18, "aa_type": "polar", "aa_type_no": 1},
    "TYR": {"aa_number": 19, "aa_type": "polar", "aa_type_no": 1},
    "VAL": {"aa_number": 20, "aa_type": "polar", "aa_type_no": 1},

    # unknown? Came up in one of the structures, haven't looked into it yet
    # TODO look into it!
    "UNK": {"aa_number": 21, "aa_type": "unkown", "aa_type_no": 1},

    # nucleic acids
    # DNA
    "DC": {"aa_number": 31, "aa_type": "unkown", "aa_type_no": 1},
    "DG": {"aa_number": 32, "aa_type": "unkown", "aa_type_no": 1},
    "DA": {"aa_number": 33, "aa_type": "unkown", "aa_type_no": 1},
    "DT": {"aa_number": 34, "aa_type": "unkown", "aa_type_no": 1},

    # RNA
    "C": {"aa_number": 41, "aa_type": "unkown", "aa_type_no": 1},
    "G": {"aa_number": 42, "aa_type": "unkown", "aa_type_no": 1},
    "A": {"aa_number": 43, "aa_type": "unkown", "aa_type_no": 1},
    "U": {"aa_number": 44, "aa_type": "unkown", "aa_type_no": 1}

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
# build_assembly = True
# build_assembly_id = 1

pdb_id = pdb_code
one_nanometre_size_in_metres = nanometre_scale * 0.1

# download the required model
if (fetch_pdb):
    pdb = atomium.fetch(pdb_id)
    output_name = pdb_id

else:
    pdb_id = molecule_name
    pdb = atomium.open(pdb_path)

    # set the molecule name for local import
    if molecule_name == "":
        molecule_name = "new_molecule"

    output_name = molecule_name

pdb_backup = pdb

assemblies = pdb.assemblies

# If true, the biological assembly will be built first and then imported.
# This can likely be translated to a set of geometry nodes to reduce computation
# time and also make it dynamic, but that will involve coding the creation of a bunch
# of nodes which will be a nightmare. Colouring all of the atoms will be problematic as
# well until the named attributes are properly released in blender 3.2, for now this
# will be fine as viewport performance for the point clouds is very good. Will increase
# the load times though depending on the assembly, as it will drastically increase the
# number of atoms to iterate over.
if build_assembly:
    # try to build the array with the given id, but if it fails (such as wrong id) continue on with the original mode.
    # TODO add a warning that will be displayed upon failure to build array.
    try:
        pdb = pdb.generate_assembly(build_assembly_id)
    except:
        print(Warning("Failed to generate biological assembly of id" +
              str(build_assembly_id)))


#pdb = atomium.open("C:\\Users\\BradyJohnston\\Desktop\\atp-frames.pdb")

# check the number of models in the file

if build_assembly:
    first_model = pdb
    n_models = 1
else:
    first_model = pdb.models[0]
    n_models = len(pdb.models)

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


def try_lookup(dict, key, value_on_fail=0):
    """
    Try looking up the value from the key in the dictionary, 
    otherwise return the value on fail so that things can keep moving
    """
    try:
        return dict[key]
    except:
        return value_on_fail


def try_append(list, value, value_on_fail=0):
    """
    Tries to append the value to the list, and adds instead the value on fail
    instead if a lookup into one of the dictionaries has failed.
    """
    try:
        list.append(value)
    except:
        list.append(value_on_fail)


def element_from_atom_name(atom_name):
    return re.findall("^[A-Z][a-z]?", atom_name)[0]


def get_element(atom):

    # try:
    #     element = atom.element
    # except:
    #     element = element_from_atom_name(atom.name)
    # if element == "X" or element == "" or element == None:
    element = element_from_atom_name(atom.name)
    return element


def get_element_num(element):
    try:
        element_number = element_dict.get(element).get("atomic_number")
    except:
        element_number = 0

    return element_number


def get_chain_char(atom):
    try:
        return atom.chain
    except:
        return "X"

# def get_aa_sequence_number(atom):
#     try:
#         atom.


for chain in first_model.chains():
    current_chain = chain.id
    for res in chain.residues():
        current_aa_id_char = res.name
        # the numbers at the end of the AA identifier "ASP.19" etc
        current_aa_sequence_number = int(
            re.findall(r"\d+", res.id.split(".")[1])[0])

        for atom in res.atoms():
            try_append(atom_id, atom.id)
            try_append(atom_location, atom.location)
            try_append(atom_element_char, get_element(atom))
            try_append(atom_element_num, get_element_num(get_element(atom)))
            try_append(atom_name_char, atom.name)
            try_append(atom_chain_char, current_chain)
            try_append(atom_aa_sequence_number, current_aa_sequence_number)
            try_append(atom_aa_id_char, current_aa_id_char)
            try_append(atom_aa_id_number, try_lookup(
                try_lookup(AA_dict, current_aa_id_char), "aa_number"))
            # try_append(atom_aa_id_number, AA_dict[current_aa_id_char]["aa_number"])
            try_append(atom_b_factor, atom.bvalue)
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

atom_chain_num = list(
    map(lambda x: int(np.where(x == unique_chains)[0]), atom_chain_char))
atom_chain_num = np.array(atom_chain_num)


unique_atoms = np.array(list(set(atom_name_char)))
unique_atoms_inds = unique_atoms.argsort()
unique_atoms = unique_atoms[unique_atoms_inds]


atom_name_num = list(
    map(lambda x: int(try_lookup(atom_name_dict, x)), atom_name_char))
atom_name_num = np.array(atom_name_num)


def create_model(name, collection, locations, bonds=[], faces=[]):
    """
    Creates a mesh with the given name in the given collection, from the supplied
    values for the locationso of vertices, and if supplied, bonds and faces.
    """
    # create a new mesh
    atom_mesh = bpy.data.meshes.new(name)
    atom_mesh.from_pydata(locations, bonds, faces)
    new_object = bpy.data.objects.new(name, atom_mesh)
    collection.objects.link(new_object)
    return new_object


# def lists_to_vec3_list(xvalue, yvalue, zvalue):
#     list_of_vectors = list(map(
#             lambda x: [xvalue[x], yvalue[x], zvalue[x]], range(len())
#         )
#     )

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
        locations = list(map(lambda x: [get_value(prop_x, x), get_value(prop_y, x), get_value(prop_z, x)], range(len(first_model.atoms())))))


def get_bond_list(model, connect_cutoff=0.35, search_distance=2):
    """
    For all atoms in the model, search for the nearby atoms given the current 
    distance cutoff, then calculate whether or not they will be bonded to their 
    nearby atoms.

    Returns a list of lists, each with two integers in them, specifying the 
    indices of the two atoms that are to be bonded.
    """

    mod = model
    mod.optimise_distances()
    for chain in mod.chains():
        for atom in chain.atoms():
            primary_radius = radii_dict[atom.element]
            nearby_atoms = atom.nearby_atoms(search_distance)
            if atom.element == "H":
                connect_adjust = -0.2
            else:
                connect_adjust = 0

            for atom2 in nearby_atoms:
                same_chain =  (atom.chain.name == atom2.chain.name)
                disulfide = ((atom.name == atom2.name) and (atom.name == 'SG'))

                if not same_chain and not disulfide:
                    continue

                # if both atoms are the sulfurs in cysteins, then treat as a 
                # disulfide bond which can be bonded from a longer distance

                if disulfide:
                    connect_adjust == 0.2
                
                secondary_radius = radii_dict[atom2.element]
                distance = atom.distance_to(atom2)
                if distance <= ((connect_cutoff + connect_adjust) + (primary_radius + secondary_radius) / 2):
                    atom.bond(atom2)

    all_atoms = mod.atoms()
    all_ids = np.array(list(map(lambda x: x.id, all_atoms)))
    inds = all_ids.argsort()
    all_ids = all_ids[inds]

    bond_list = []

    for atom in all_atoms:
        for atom2 in atom.bonded_atoms:
            bond_list.append([int(np.where(atom.id == all_ids)[0]), int(
                np.where(atom2.id == all_ids)[0])])

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


def get_frame_bvalue(frame):
    """
    Returns a numpy array of all of the atom bvalue from the given frame. 
    Importantly it orders them according to their atom numbering to sync the frames.
    """
    all_atoms = frame.atoms()
    atom_id = list(map(lambda x: x.id, all_atoms))
    atom_bvalue = list(map(lambda x: x.bvalue, all_atoms))

    atom_id = np.array(atom_id)
    inds = atom_id.argsort()
    atom_id = atom_id[inds]
    atom_bvalue = np.array(atom_bvalue)
    atom_bvalue = atom_bvalue[inds]

    return atom_bvalue


def get_model_element_number(model):
    """
    Returns a numpy array of all of the atom bvalue from the given frame. 
    Importantly it orders them according to their atom numbering to sync the frames.
    """
    def try_element_number(element):
        try:
            return element_dict[element]["atomic_number"]
        except:
            return 3

    all_atoms = model.atoms()
    atom_id = list(map(lambda x: x.id, all_atoms))
    atom_element = list(map(lambda x: get_element(x), all_atoms))

    atom_element_number = list(
        map(lambda x: try_element_number(x), atom_element))

    atom_id = np.array(atom_id)
    inds = atom_id.argsort()
    atom_id = atom_id[inds]
    atom_element_number = np.array(atom_element_number)
    atom_element_number = atom_element_number[inds]

    return atom_element_number


def get_model_is_sidechain(model):
    """
    Returns a numpy array of all of the atom bvalue from the given frame. 
    Importantly it orders them according to their atom numbering to sync the frames.
    """
    def try_is_sidechain(atom):
        try:
            return int(atom.is_side_chain)
        except:
            return 0

    all_atoms = model.atoms()
    atom_id = list(map(lambda x: x.id, all_atoms))
    atom_is_sidechain = list(map(lambda x: try_is_sidechain(x), all_atoms))

    atom_id = np.array(atom_id)
    inds = atom_id.argsort()
    atom_id = atom_id[inds]
    atom_is_sidechain = np.array(atom_is_sidechain)
    atom_is_sidechain = atom_is_sidechain[inds]

    return atom_is_sidechain


def get_model_is_backbone(model):
    """
    Returns a numpy array of all of the atom bvalue from the given frame. 
    Importantly it orders them according to their atom numbering to sync the frames.
    """
    def try_is_backbone(atom):
        try:
            return int(atom.is_backbone)
        except:
            return 0

    all_atoms = model.atoms()
    atom_id = list(map(lambda x: x.id, all_atoms))
    atom_is_backbone = list(map(lambda x: try_is_backbone(x), all_atoms))

    atom_id = np.array(atom_id)
    inds = atom_id.argsort()
    atom_id = atom_id[inds]
    atom_is_backbone = np.array(atom_is_backbone)
    atom_is_backbone = atom_is_backbone[inds]

    return atom_is_backbone


def get_model_is_ca(model):
    """
    Returns numpy array of True / False if atom name is equal to CA.
    """

    all_atoms = model.atoms()
    try:
        atom_id = list(map(lambda x: x.id, all_atoms))
        atom_is_CA = list(map(lambda x: x.name == "CA", all_atoms))

        atom_id = np.array(atom_id)
        inds = atom_id.argsort()
        atom_is_CA = np.array(atom_is_CA)[inds]
    
    except:
        atom_is_CA = []
        for i in len(all_atoms):
            atom_is_CA.append(False)

    return atom_is_CA

# See if there is a collection called "Molecular Nodes", if so, set it to be the parent
# collection, otherwise create one and link it to the scene collection.
try:
    parent_coll = bpy.data.collections['MolecularNodes']
    parent_coll.name == "MolecularNodes"
    # make the collection active, for creating and disabling new
    bpy.context.view_layer.active_layer_collection = bpy.context.view_layer.layer_collection.children[
        'MolecularNodes']
except:
    parent_coll = bpy.data.collections.new('MolecularNodes')
    bpy.context.scene.collection.children.link(parent_coll)
    # make the collection active, for creating and disabling new
    bpy.context.view_layer.active_layer_collection = bpy.context.view_layer.layer_collection.children[
        'MolecularNodes']


# create new collection that will house the data, link it to the parent collection
col = bpy.data.collections.new(output_name)
parent_coll.children.link(col)

col_properties = bpy.data.collections.new(output_name + "_properties")
col.children.link(col_properties)

# If create_bonds selected, generate a list of vertex pairs that will be the bonds for the atomic mesh,
# else return an empty list that will make no edges when passed to create_model()
if create_bonds:
    bonds = get_bond_list(first_model, connect_cutoff=connect_cutoff)
else:
    bonds = []

# create the first model, that will be the actual atomic model the user will interact with and display
base_model = create_model(
    name=output_name,
    collection=col,
    locations=get_frame_positions(first_model) * one_nanometre_size_in_metres,
    bonds=bonds
)

# Creat the different models that will encode the various properties into
# the XYZ locations of ther vertices.
create_properties_model(
    name=output_name + "_properties_1",
    collection=col_properties,
    prop_x=get_model_element_number(first_model),
    prop_y=atom_chain_num + 1,  # to have the first chain be indexed from 1
    prop_z=atom_name_num
)
create_properties_model(
    name=output_name + "_properties_2",
    collection=col_properties,
    prop_x=atom_aa_sequence_number,
    prop_y=atom_id,
    prop_z=atom_aa_id_number
)

create_properties_model(
    name=output_name + "_properties_3",
    collection=col_properties,
    prop_x=get_frame_bvalue(first_model),
    prop_y=get_model_is_backbone(first_model),
    prop_z=get_model_is_ca(first_model)
)

# hide the created properties collection
bpy.context.layer_collection.children[col.name].children[col_properties.name].exclude = True


# create the frames
if (n_models > 1):
    col_frames = bpy.data.collections.new(output_name + "_frames")
    col.children.link(col_frames)
    # for each model in the pdb, create a new object and add it to the frames collection
    # testing out the addition of points that represent the bfactors. You can then in theory
    # use the modulo of the index to be able to pick either the position or the bvalue for
    # each frame in the frames collection.
    for frame in pdb_backup.models:
        atom_locations = get_frame_positions(
            frame) * one_nanometre_size_in_metres
        # need to make sure that is a list of 3-element vectors to encode to the XYZ positions
        # of the vertices
        atom_bvalue = list(map(lambda x: [x, 0, 0], get_frame_bvalue(frame)))
        model_locations = list(atom_locations) + atom_bvalue
        create_model(
            name="frame_" + output_name,
            collection=col_frames,
            locations=model_locations
        )

    # hide the created frames collection
    bpy.context.layer_collection.children[col.name].children[col_frames.name].exclude = True
