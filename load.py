import bpy
import numpy as np
# import biotite.structure as struc
# import biotite.structure.io.pdb as pdb
# import biotite.structure.io.pdbx as pdbx
# import biotite.structure.io.mmtf as mmtf
# import biotite.database.rcsb as rcsb
from .tools import coll_mn
import warnings
from . import data

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
        # warnings.warn("Unable to create attribute: " + name, bpy.)
        return None

def create_molecule(mol_array, mol_name, center_molecule = False, del_solvent = False, include_bonds = True):
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

    if include_bonds:
        bonds = mol_array.bonds.as_array()
        mol_object = create_object(name = mol_name, collection = coll_mn(), locations = locations, bonds = bonds[:, [0,1]])
    else:
        mol_object = create_object(name = mol_name, collection = coll_mn(), locations = locations)

    # compute the attributes as numpy arrays for the addition of them to the points of the structure
    atomic_number = np.fromiter(map(lambda x: data.elements.get(x, {'atomic_number': 0}).get("atomic_number"), np.char.title(mol_array.element)), dtype = np.int)
    res_id = mol_array.res_id
    res_name = np.fromiter(map(lambda x: data.amino_acids.get(x, {'aa_number': 0}).get('aa_number'), np.char.upper(mol_array.res_name)), dtype = np.int)
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
    

    # add bonds attribute to the model if there was bond information
    if include_bonds:
        try:
            bond_type = bonds[:, 2].copy(order = 'C') # the .copy(order = 'C') is to fix a weird ordering issue with the resulting array
            add_attribute(mol_object, 'bond_type', bond_type,"INT", "EDGE")
        except:
            warnings.warn('Unable to add bond information to the molecule.')

    
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
        coll_mn().children.link(coll_frames)
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
        bpy.context.view_layer.layer_collection.children[coll_mn().name].children[coll_frames.name].exclude = True
    
    # add custom properties for the object, such as number of chains, biological assemblies etc
    try:
        mol_object['chain_id_unique'] = list(np.unique(mol_array.chain_id))
    except:
        warnings.warn('No chain information detected.')
    
    return mol_object
