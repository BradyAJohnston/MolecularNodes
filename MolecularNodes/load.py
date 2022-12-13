import bpy
import numpy as np
from .tools import coll_mn
import warnings
from . import data
from . import assembly
from . import nodes

def molecule_rcsb(pdb_code, 
                  center_molecule=False, 
                  del_solvent=True, 
                  include_bonds=True, 
                  starting_style=0, 
                  setup_nodes=True
                  ):
    
    mol, file = open_structure_rcsb(pdb_code = pdb_code, include_bonds=include_bonds)
    mol_object, coll_frames = create_molecule(
        mol_array = mol,
        mol_name = pdb_code,
        center_molecule = center_molecule,
        del_solvent = del_solvent, 
        include_bonds = include_bonds
        )
    
    if setup_nodes:
        nodes.create_starting_node_tree(
            obj = mol_object, 
            coll_frames=coll_frames, 
            starting_style = starting_style
            )
    
    mol_object['bio_transform_dict'] = file['bioAssemblyList']
    
    return mol_object

def molecule_local(file_path, 
                   mol_name="Name",
                   include_bonds=True, 
                   center_molecule=True, 
                   del_solvent=True, 
                   default_style=0, 
                   setup_nodes=True
                   ): 
    import biotite.structure as struc
    
    
    import os
    file_path = os.path.abspath(file_path)
    file_ext = os.path.splitext(file_path)[1]
    
    if file_ext == '.pdb':
        mol, file = open_structure_local_pdb(file_path, include_bonds)
        transforms = assembly.get_transformations_pdb(file)
    elif file_ext == '.pdbx' or file_ext == '.cif':
        mol, file = open_structure_local_pdbx(file_path, include_bonds)
        try:
            transforms = assembly.get_transformations_pdbx(file)
        except:
            transforms = None
            # self.report({"WARNING"}, message='Unable to parse biological assembly information.')
    
    # if include_bonds chosen but no bonds currently exist (mol.bonds is None)
    # then attempt to find bonds by distance
    if include_bonds and not mol.bonds:
        mol.bonds = struc.connect_via_distances(mol[0], inter_residue=True)
    
    mol_object, coll_frames = create_molecule(
        mol_array = mol,
        mol_name = mol_name,
        center_molecule = center_molecule,
        del_solvent = del_solvent, 
        include_bonds = include_bonds
        )
    
        
    
    # setup the required initial node tree on the object 
    if setup_nodes:
        nodes.create_starting_node_tree(
            obj = mol_object,
            coll_frames = coll_frames,
            starting_style = default_style
            )
    
    # if transforms:
        # mol_object['bio_transform_dict'] = (transforms)
        # mol_object['bio_transnform_dict'] = 'testing'
        
    return mol_object


def open_structure_rcsb(pdb_code, include_bonds = True):
    import biotite.structure.io.mmtf as mmtf
    import biotite.database.rcsb as rcsb
    
    file = mmtf.MMTFFile.read(rcsb.fetch(pdb_code, "mmtf"))
    
    # returns a numpy array stack, where each array in the stack is a model in the 
    # the file. The stack will be of length = 1 if there is only one model in the file
    mol = mmtf.get_structure(file, extra_fields = ["b_factor", "charge"], include_bonds = include_bonds) 
    return mol, file


def open_structure_local_pdb(file_path, include_bonds = True):
    import biotite.structure.io.pdb as pdb
    
    file = pdb.PDBFile.read(file_path)
    
    # returns a numpy array stack, where each array in the stack is a model in the 
    # the file. The stack will be of length = 1 if there is only one model in the file
    mol = pdb.get_structure(file, extra_fields = ['b_factor', 'charge'], include_bonds = include_bonds)
    return mol, file

def open_structure_local_pdbx(file_path, include_bonds = True):
    import biotite.structure as struc
    import biotite.structure.io.pdbx as pdbx
    
    file = pdbx.PDBxFile.read(file_path)
    
    # returns a numpy array stack, where each array in the stack is a model in the 
    # the file. The stack will be of length = 1 if there is only one model in the file
    mol  = pdbx.get_structure(file, extra_fields = ['b_factor', 'charge'])
    # pdbx doesn't include bond information apparently, so manually create
    # them here if requested
    if include_bonds:
        mol[0].bonds = struc.bonds.connect_via_residue_names(mol[0], inter_residue = True)
    return mol, file

def create_object(name, collection, locations, bonds=[]):
    """
    Creates a mesh with the given name in the given collection, from the supplied
    values for the locations of vertices, and if supplied, bonds as edges.
    """
    # create a new mesh
    mol_mesh = bpy.data.meshes.new(name)
    mol_mesh.from_pydata(locations, bonds, faces=[])
    mol_object = bpy.data.objects.new(name, mol_mesh)
    collection.objects.link(mol_object)
    return mol_object

def add_attribute(object, name, data, type = "FLOAT", domain = "POINT", add = True):
    if not add:
        return None
    try:
        attribute = object.data.attributes.new(name, type, domain)
        attribute.data.foreach_set('value', data)
        return True
    except:
        warnings.warn("Unable to create attribute: " + name)
        return None

def create_molecule(mol_array, mol_name, center_molecule = False, del_solvent = False, include_bonds = False, collection = None):
    import biotite.structure as struc
    
    if np.shape(mol_array)[0] > 1:
        mol_frames = mol_array
    else:
        mol_frames = None
    
    mol_array = mol_array[0]
    
    # remove the solvent from the structure if requested
    if del_solvent:
        mol_array = mol_array[np.invert(struc.filter_solvent(mol_array))]

    world_scale = 0.01
    locations = mol_array.coord * world_scale
    
    centroid = np.array([0, 0, 0])
    if center_molecule:
        centroid = struc.centroid(mol_array) * world_scale
    

    # subtract the centroid from all of the positions to localise the molecule on the world origin
    if center_molecule:
        locations = locations - centroid

    if not collection:
        collection = coll_mn()
    
    if include_bonds:
        bonds = mol_array.bonds.as_array()
        mol_object = create_object(name = mol_name, collection = collection, locations = locations, bonds = bonds[:, [0,1]])
    else:
        mol_object = create_object(name = mol_name, collection = collection, locations = locations)

    # compute the attributes as numpy arrays for the addition of them to the points of the structure
    # TODO find a way to do this nicer, and with more control when something fails
    atomic_number = np.fromiter(map(lambda x: data.elements.get(x, {'atomic_number': -1}).get("atomic_number"), np.char.title(mol_array.element)), dtype = np.int)
    res_id = mol_array.res_id
    res_name = np.fromiter(map(lambda x: data.residues.get(x, {'res_name_num': -1}).get('res_name_num'), np.char.upper(mol_array.res_name)), dtype = np.int)
    chain_id = np.searchsorted(np.unique(mol_array.chain_id), mol_array.chain_id)
    vdw_radii =  np.fromiter(map(struc.info.vdw_radius_single, mol_array.element), dtype=np.float) * world_scale
    is_alpha = np.fromiter(map(lambda x: x == "CA", mol_array.atom_name), dtype = np.bool)
    is_solvent = struc.filter_solvent(mol_array)
    is_backbone = (struc.filter_backbone(mol_array) | np.isin(mol_array.atom_name, ["P", "O5'", "C5'", "C4'", "C3'", "O3'"]))
    is_nucleic = struc.filter_nucleotides(mol_array)
    is_peptide = struc.filter_canonical_amino_acids(mol_array)
    is_hetero = mol_array.hetero
    is_carb = struc.filter_carbohydrates(mol_array)
    b_factor = mol_array.b_factor

    # Add information about the bond types to the model on the edge domain
    # Bond types: 'ANY' = 0, 'SINGLE' = 1, 'DOUBLE' = 2, 'TRIPLE' = 3, 'QUADRUPLE' = 4
    # 'AROMATIC_SINGLE' = 5, 'AROMATIC_DOUBLE' = 6, 'AROMATIC_TRIPLE' = 7
    # https://www.biotite-python.org/apidoc/biotite.structure.BondType.html#biotite.structure.BondType
    if include_bonds:
        try:
            add_attribute(
                object = mol_object, 
                name = 'bond_type', 
                data = bonds[:, 2].copy(order = 'C'), # the .copy(order = 'C') is to fix a weird ordering issue with the resulting array
                type = "INT", 
                domain = "EDGE"
                )
        except:
            warnings.warn('Unable to add bond types to the molecule.')

    
    # these are all of the attributes that will be added to the structure
    # TODO add capcity for selection of particular attributes to include / not include to potentially
    # boost performance, unsure if actually a good idea of not. Need to do some testing.
    attributes = (
        {'name': 'res_id',          'value': res_id,              'type': 'INT',     'domain': 'POINT'},
        {'name': 'res_name',        'value': res_name,            'type': 'INT',     'domain': 'POINT'},
        {'name': 'atomic_number',   'value': atomic_number,       'type': 'INT',     'domain': 'POINT'},
        {'name': 'b_factor',        'value': b_factor,            'type': 'FLOAT',   'domain': 'POINT'},
        {'name': 'vdw_radii',       'value': vdw_radii,           'type': 'FLOAT',   'domain': 'POINT'},
        {'name': 'chain_id',        'value': chain_id,            'type': 'INT',     'domain': 'POINT'},
        {'name': 'is_backbone',     'value': is_backbone,         'type': 'BOOLEAN', 'domain': 'POINT'},
        {'name': 'is_alpha_carbon', 'value': is_alpha,            'type': 'BOOLEAN', 'domain': 'POINT'},
        {'name': 'is_solvent',      'value': is_solvent,          'type': 'BOOLEAN', 'domain': 'POINT'},
        {'name': 'is_nucleic',      'value': is_nucleic,          'type': 'BOOLEAN', 'domain': 'POINT'},
        {'name': 'is_peptide',      'value': is_peptide,          'type': 'BOOLEAN', 'domain': 'POINT'},
        {'name': 'is_hetero',       'value': is_hetero,           'type': 'BOOLEAN', 'domain': 'POINT'},
        {'name': 'is_carb',         'value': is_carb,             'type': 'BOOLEAN', 'domain': 'POINT'}
    )
    
    # assign the attributes to the object
    for att in attributes:
        try:
            add_attribute(mol_object, att['name'], att['value'], att['type'], att['domain'])
        except:
            warnings.warn('Unable to add ' + att['name'] + ' to the ' + att['domain'] + ' domain.')
    
    if mol_frames:
        # create the frames of the trajectory in their own collection to be disabled
        coll_frames = bpy.data.collections.new(mol_object.name + "_frames")
        collection.children.link(coll_frames)
        counter = 0
        for frame in mol_frames:
            obj_frame = create_object(
                name = mol_object.name + '_frame_' + str(counter), 
                collection=coll_frames, 
                locations= frame.coord * world_scale - centroid
            )
            try:
                add_attribute(obj_frame, 'b_factor', frame.b_factor)
            except:
                pass
            counter += 1
        
        # disable the frames collection so it is not seen
        bpy.context.view_layer.layer_collection.children[collection.name].children[coll_frames.name].exclude = True
    else:
        coll_frames = None
    
    # add custom properties to the actual blender object, such as number of chains, biological assemblies etc
    # currently biological assemblies can be problematic to holding off on doing that
    try:
        mol_object['chain_id_unique'] = list(np.unique(mol_array.chain_id))
    except:
        warnings.warn('No chain information detected.')
    
    return mol_object, coll_frames
