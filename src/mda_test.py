import bpy
from logging import warning
from gc import collect
import sys
import os
import site
from time import time
import numpy as np


def verify_user_sitepackages(mda_path):
    usersitepackagespath = site.getusersitepackages()

    if os.path.exists(usersitepackagespath) and usersitepackagespath not in sys.path:
        sys.path.append(usersitepackagespath)
    if os.path.exists(mda_path) and mda_path not in sys.path:
        sys.path.append(mda_path)


verify_user_sitepackages(mdanalysis_dir_location)

try:
    import MDAnalysis as mda
except:
    warning("Unable to Import MDAnalysis")

dict_elements = {
    "H"  :  1,
    "He" :  2,
    "Li" :  3,
    "Be" :  4,
    "B"  :  5,
    "C"  :  6,
    "N"  :  7,
    "O"  :  8,
    "F"  :  9,
    "Ne" : 10,
    "Na" : 11,
    "Mg" : 12,
    "Al" : 13,
    "Si" : 14,
    "P"  : 15,
    "S"  : 16,
    "Cl" : 17,
    "Ar" : 18,
    "K"  : 19,
    "Ca" : 20,
    "Sc" : 21,
    "Ti" : 22,
    "V"  : 23,
    "Cr" : 24,
    "Mn" : 25,
    "Fe" : 26,
    "Co" : 27,
    "Ni" : 28,
    "Cu" : 29,
    "Zn" : 30,
    "Ga" : 31,
    "Ge" : 32,
    "As" : 33,
    "Se" : 34,
    "Br" : 35,
    "Kr" : 36,
    "Rb" : 37,
    "Sr" : 38,
    "Y"  : 39,
    "Zr" : 40,
    "Nb" : 41,
    "Mo" : 42,
    "Tc" : 43,
    "Ru" : 44,
    "Rh" : 45,
    "Pd" : 46,
    "Ag" : 47,
    "Cd" : 48,
    "In" : 49,
    "Sn" : 50,
    "Sb" : 51,
    "Te" : 52,
    "I"  : 53,
    "Xe" : 54,
    "Cs" : 55,
    "Ba" : 56,
    "La" : 57,
    "Ce" : 58,
    "Pr" : 59,
    "Nd" : 60,
    "Pm" : 61,
    "Sm" : 62,
    "Eu" : 63,
    "Gd" : 64,
    "Tb" : 65,
    "Dy" : 66,
    "Ho" : 67,
    "Er" : 68,
    "Tm" : 69,
    "Yb" : 70,
    "Lu" : 71,
    "Hf" : 72,
    "Ta" : 73,
    "W"  : 74,
    "Re" : 75,
    "Os" : 76,
    "Ir" : 77,
    "Pt" : 78,
    "Au" : 79,
    "Hg" : 80,
    "Tl" : 81,
    "Pb" : 82,
    "Bi" : 83,
    "Po" : 84,
    "At" : 85,
    "Rn" : 86,
    "Fr" : 87,
    "Ra" : 88,
    "Ac" : 89,
    "Th" : 90,
    "Pa" : 91,
    "U"  : 92,
    "Np" : 93,
    "Pu" : 94,
    "Am" : 95,
    "Cm" : 96,
    "Bk" : 97,
    "Cf" : 98,
    "Es" : 99,
    "Fm" : 100,
    "Md" : 101,
    "No" : 102,
    "Lr" : 103,
    "Vac": 104
}


# See if there is a collection called "Molecular Nodes", if so, set it to be the parent
# collection, otherwise create one and link it to the scene collection.

#top = "G:\\bjohnston\\BAJ_01\\box.gro"
#
# traj = "G:\\bjohnston\\BAJ_01\\md_nojump.xtc"
# file_top = "C:\\Users\\BradyJohnston\\Downloads\\6vsb_2_2_2_traj_xtc\\last_frame_nos.pdb"
# file_traj = "C:\\Users\\BradyJohnston\\Downloads\\6vsb_2_2_2_traj_xtc\\trj_nos.xtc"

# setup the universe using MDAnalysis, for use later in the script
u = mda.Universe(file_top, file_traj)

# output_name = "spike"
output_name = molecule_name

# try to get the collection, returns None if collection doesn't exist
parent_coll = bpy.data.collections.get('MolecularNodes')

# if the collection doesn't exist, create it
if not parent_coll:
    parent_coll = bpy.data.collections.new('MolecularNodes')
    bpy.context.scene.collection.children.link(parent_coll)

# make the MolecularNodes collection active
bpy.context.view_layer.active_layer_collection = bpy.context.view_layer.layer_collection.children['MolecularNodes']

# try:
#     parent_coll = bpy.data.collections['MolecularNodes']
#     parent_coll.name == "MolecularNodes"
#     # make the collection active, for creating and disabling new
#     bpy.context.view_layer.active_layer_collection = bpy.context.view_layer.layer_collection.children[
#         'MolecularNodes']
# except:
#     parent_coll = bpy.data.collections.new('MolecularNodes')
#     bpy.context.scene.collection.children.link(parent_coll)
#     # make the collection active, for creating and disabling new
#     bpy.context.view_layer.active_layer_collection = bpy.context.view_layer.layer_collection.children[
#         'MolecularNodes']


# create new collection that will house the data, link it to the parent collection
col = bpy.data.collections.new(output_name)
parent_coll.children.link(col)

# create the properties collection
col_properties = bpy.data.collections.new(output_name + "_properties")
col.children.link(col_properties)


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


def create_properties_model(name, collection, prop_x, prop_y, prop_z, n_atoms):
    """
    Creates a mesh that will act as a look up table for properties about the atoms
    in the actual mesh that will be created.
    """
    if n_atoms == None:
        n_atoms = len(prop_x)

    def get_value(vec, x):
        try:
            return vec[x]
        except:
            return 0

    list_list = list(map(
        lambda x: [
            get_value(prop_x, x),
            get_value(prop_y, x),
            get_value(prop_z, x)
        ], range(n_atoms - 1)
    ))

    create_model(
        name=name,
        collection=collection,
        locations=list_list
    )


def create_prop_AtomicNum_ChainNum_NameNum(univ, name, collection):
    n = univ.atoms.n_atoms

    # get the atomic numbers for the atoms
    try:
        prop_elem_num = list(
            map(lambda x: dict_elements.get(x, 0), u.atoms.elements))

    except:
        try:
            prop_elem = list(
                map(lambda x: mda.topology.guessers.guess_atom_element(x), univ.atoms.names))
            prop_elem_num = list(
                map(lambda x: dict_elements.get(x, 0), prop_elem))
        except:
            prop_elem = []
            for i in range(n - 1):
                prop_elem.append(0)

    # get the chain numbers for the atoms
    try:
        prop_chain = univ.atoms.chainIDs
        unique_chains = np.array(list(set(prop_chain)))
        prop_chain_num = np.array(
            list(map(lambda x: np.where(x == unique_chains)[0][0], prop_chain)))
    except:
        unique_chains = ["No Chains Found"]
        prop_chain_num = []
        for i in range(n - 1):
            prop_chain_num.append(0)

    # get the name numbers
    try:
        prop_name = univ.atoms.names
        unique_names = np.array(list(set(prop_name)))
        prop_name_num = np.array(
            list(map(lambda x: np.where(x == unique_names)[0][0], prop_name)))
    except:
        prop_name = []
        for i in range(n - 1):
            prop_name.append(0)

    create_properties_model(
        name=name,
        collection=collection,
        prop_x=prop_elem_num,
        prop_y=prop_chain_num,
        prop_z=prop_name_num,
        n_atoms=n
    )

    return list(unique_chains)


def create_prop_aaSeqNum_atomNum_atomAAIDNum(univ, name, collection):
    n = univ.atoms.n_atoms

    # get the AA sequence numebrs
    try:
        prop_aaSeqNum = univ.atoms.resnums
    except:
        prop_aaSeqNum = []
        for i in range(n - 1):
            prop_aaSeqNum.append(0)

    # get the atom indices
    try:
        prop_atomNum = univ.atoms.ids
    except:
        prop_atomNum = range(1, n + 1)

    # get the residue names (AA names etc)
    try:
        resnames = univ.atoms.resnames
        unique_resnames = np.array(list(set(resnames)))
        prop_aa_ID_num = list(map(
            lambda x: np.where(x == unique_resnames)[0][0], resnames
        ))

    except:
        prop_aa_ID_num = []
        for i in range(n - 1):
            prop_aa_ID_num.append(0)

    create_properties_model(
        name=name,
        collection=collection,
        prop_x=prop_aaSeqNum,
        prop_y=prop_atomNum,
        prop_z=prop_aa_ID_num,
        n_atoms=n
    )


def create_prop_bvalue_isBackbone_isCA(univ, name, collection):
    n = univ.atoms.n_atoms

    # setup bvalue properties, likely that most simulations won't actually
    # have bfactors, but here anyway to capture them if they do.
    try:
        prop_bvalue = univ.atoms.tempfactors
    except:
        prop_bvalue = []
        for i in range(n - 1):
            prop_bvalue.append(0)

    # setup isBackbone properties, selects backbone for nucleic and protein
    try:
        prop_is_backbone = np.isin(univ.atoms.ix, univ.select_atoms(
            'backbone or nucleicbackbone').ix.astype(int))
    except:
        prop_is_backbone = []
        for i in range(n - 1):
            prop_bvalue.append(0)

    try:
        # compare the indices against a subset of indices for only the alpah carbons,
        # convert it to ingeger of 0 = False and 1 = True
        prop_is_CA = np.isin(
            univ.atoms.ix, univ.select_atoms("name CA").ix).astype(int)
    except:
        prop_is_CA = []
        for i in range(n - 1):
            prop_is_CA.append(0)

    create_properties_model(
        name=name,
        collection=collection,
        prop_x=prop_bvalue,
        prop_y=prop_is_backbone,
        prop_z=prop_is_CA,
        n_atoms=n
    )


# create the first model, that will be the actual atomic model the user will interact with and display
base_model = create_model(
    name=output_name,
    collection=col,
    locations=u.atoms.positions * 0.1,
    bonds=[]
)


# create the models that will hold the properties associated with each atom
unique_names = create_prop_AtomicNum_ChainNum_NameNum(
    univ=u,
    name=output_name + "_properties_1",
    collection=col_properties
)
create_prop_aaSeqNum_atomNum_atomAAIDNum(
    univ=u,
    name=output_name + "_properties_2",
    collection=col_properties
)

create_prop_bvalue_isBackbone_isCA(
    univ=u,
    name=output_name + "_properties_3",
    collection=col_properties
)

# create the frames

col_frames = bpy.data.collections.new(output_name + "_frames")
col.children.link(col_frames)
# for each model in the pdb, create a new object and add it to the frames collection
# testing out the addition of points that represent the bfactors. You can then in theory
# use the modulo of the index to be able to pick either the position or the bvalue for
# each frame in the frames collection.


def create_frames(universe, collection, start=1, end=50000, time_step=100, name=output_name, nm_scale=0.1):
    """
    From the given universe, add frames to the given collection from the start till the end, along the given time
    step for each frame.
    """
    counter = 1
    for ts in universe.trajectory:
        if counter % time_step == 0 and counter > start and counter < end:
            create_model(
                name=output_name + "_frame_" + str(counter),
                collection=collection,
                locations=universe.atoms.positions * nm_scale
            )
        counter += 1


# create the frames from the given universe, only along the given timesteps
create_frames(
    universe=u,
    collection=col_frames,
    start=md_frame_start,
    time_step=md_frame_interval,
    end=md_frame_end,
    name=output_name,
    nm_scale=0.1
)

# hide the created frames collection and the properties collection
bpy.context.layer_collection.children[col.name].children[col_frames.name].exclude = True
bpy.context.layer_collection.children[col.name].children[col_properties.name].exclude = True

# # try to get the Molecular Nodes modifier and select it, if not create one and select it


# def create_starting_node_tree(collection_of_properties, obj):
#     try:
#         node_mod = obj.modifiers['MolecularNodes']
#     except:
#         node_mod = None

#     if node_mod == None:
#         node_mod = obj.modifiers.new("MolecularNodes", "NODES")
#         obj.modifiers.active = node_mod
#     else:
#         obj.modifiers.active = node_mod

#     node_mod.node_group.name = "MOL_" + str(output_name)

#     node_input = node_mod.node_group.nodes['Group Input']
#     node_input.location = [-200, 0]
#     node_output = node_mod.node_group.nodes['Group Output']
#     node_output.location = [600, 0]

#     # create an empty node group and link it to the atomic properties node group
#     new_node_group = node_mod.node_group.nodes.new("GeometryNodeGroup")
#     new_node_group.node_tree = bpy.data.node_groups["MOL_atomic_properties"]
#     new_node_group.inputs['Properties'].default_value = collection_of_properties
#     new_node_group.location = [0, 0]
#     # resize the newly created node to be a bit wider
#     node_mod.node_group.nodes[-1].width = 200

#     colour_node_group = node_mod.node_group.nodes.new("GeometryNodeGroup")
#     colour_node_group.node_tree = bpy.data.node_groups["MOL_style_colour"]
#     colour_node_group.location = [300, 0]
#     node_mod.node_group.nodes[-1].width = 200

#     link = node_mod.node_group.links.new

#     link(node_input.outputs['Geometry'], new_node_group.inputs['Atoms'])
#     link(new_node_group.outputs['Atoms'], colour_node_group.inputs['Atoms'])
#     link(new_node_group.outputs['atomic_number'],
#          colour_node_group.inputs['atomic_number'])
#     link(colour_node_group.outputs['Atoms'], node_output.inputs['Geometry'])

#     node_mod.node_group.outputs.new("NodeSocketColor", "Colour")
#     link(colour_node_group.outputs['Colour'], node_output.inputs['Colour'])
#     node_mod['Output_2_attribute_name'] = "Colour"

#     mat = create_starting_material()

#     colour_node_group.inputs['Material'].default_value = mat


# def create_starting_material():
#     try:
#         mat = bpy.data.materials['MOL_atomic_material']
#         return mat
#     except:
#         mat = None

#     if mat == None:
#         mat = bpy.data.materials.new('MOL_atomic_material')

#     mat.use_nodes = True
#     node_att = mat.node_tree.nodes.new("ShaderNodeAttribute")

#     node_att.attribute_name = "Colour"
#     node_att.location = [-300, 200]

#     mat.node_tree.links.new(
#         node_att.outputs['Color'], mat.node_tree.nodes['Principled BSDF'].inputs['Base Color'])

#     return mat


# # create_starting_material(base_model)
# create_starting_node_tree(
#     collection_of_properties=col_properties, obj=base_model)

# base_model.select_set(True)
# bpy.context.view_layer.objects.active = base_model

# # bpy.context.active_object.data.materials.append(mat)
