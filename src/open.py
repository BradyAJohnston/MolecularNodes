import bpy
import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdb as pdb
import biotite.structure.io.pdbx as pdbx
import biotite.structure.io.mmtf as mmtf
import biotite.database.rcsb as rcsb
from .tools import mn_collection
import warnings

def open_structure_rcsb(pdb_code, include_bonds = True):
    file = mmtf.MMTFFile.read(rcsb.fetch(pdb_code, "mmtf"))
    mol = mmtf.get_structure(file, model = 1, extra_fields = ["b_factor", "charge"], include_bonds = include_bonds)
    return mol, file


def open_structure_local_pdb(file_path, include_bonds = True):
    file = pdb.PDBFile.read(file_path)
    mol = pdb.get_structure(file, model = 1, extra_fields = ['b_factor', 'charge'], include_bonds = include_bonds)
    return mol, file

def open_structure_local_pdbx(file_path, include_bonds = True):
    file = pdbx.PDBxFile.read(file_path)
    mol  = pdbx.get_structure(file, model = 1, extra_fields = ['b_factor', 'charge'])
    # pdbx doesn't include bond information apparently, so manually create
    # them here if requested
    if include_bonds:
        mol.bonds = struc.bonds.connect_via_residue_names(mol, inter_residue = True)
    return mol, file



def create_model(name, collection, locations, bonds=[]):
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
        return
    try:
        attribute = object.data.attributes.new(name, type, domain)
        attribute.data.foreach_set('value', data)
    except:
        # warnings.warn("Unable to create attribute: " + name, bpy.)
        return

def MOL_import_mol(mol, mol_name, center_molecule = False, del_solvent = False, include_bonds = True):
    import bpy
    import numpy as np
    import biotite.structure as struc
    import biotite.database.rcsb as rcsb
    import biotite.structure.io.mmtf as mmtf
    from . import packages
    from . import dict

    packages.install_packages()
    packages.verify()

    # remove the solvent from the structure if requested
    if del_solvent:
        mol = mol[np.invert(struc.filter_solvent(mol))]

    world_scale = 0.01
    locations = mol.coord * world_scale
    centroid = struc.centroid(mol) * world_scale

    # subtract the centroid from all of the positions to localise the molecule on the world origin
    if center_molecule:
        locations = locations - centroid

    if include_bonds:
        bonds = mol.bonds.as_array()
        mol_object = create_model(name = mol_name, collection = mn_collection(), locations = locations, bonds = bonds[:, [0,1]])
    else:
        mol_object = create_model(name = mol_name, collection = mn_collection(), locations = locations)

    # compute some of the attributes before assignment
    atomic_number = np.fromiter(map(lambda x: dict.elements.get(x, {'atomic_number': 0}).get("atomic_number"), np.char.title(mol.element)), dtype = np.int)
    res_name = np.fromiter(map(lambda x: dict.amino_acids.get(x, {'aa_number': 0}).get('aa_number'), np.char.upper(mol.res_name)), dtype = np.int)
    chain_id = np.searchsorted(np.unique(mol.chain_id), mol.chain_id)
    vdw_radii =  np.fromiter(map(struc.info.vdw_radius_single, mol.element), dtype=np.float) * world_scale
    is_alpha = np.fromiter(map(lambda x: x == "CA", mol.atom_name), dtype = np.bool)
    is_solvent = struc.filter_solvent(mol)
    is_backbone = (struc.filter_backbone(mol) | np.isin(mol.atom_name, ["P", "O5'", "C5'", "C4'", "C3'", "O3'"]))
    is_nucleic = struc.filter_nucleotides(mol)
    is_peptide = struc.filter_canonical_amino_acids(mol)
    is_carb = struc.filter_carbohydrates(mol)
    

    if include_bonds:
        bond_type = bonds[:, 2].copy(order = 'C') # the .copy(order = 'C') is to fix a weird ordering issue with the resulting array
        add_attribute(mol_object, 'bond_type', bond_type,"INT", "EDGE")

    
    # these are all of the attributes that will be added to the structure
    # TODO add capcity for selection of particular attributes to include / not include to potentially
    # boost performance, unsure if actually a good idea of not. Need to do some testing.
    attributes = (
        {'name': 'res_id',          'value': mol.res_id,    'type': 'INT',     'domain': 'POINT'},
        {'name': 'res_name',        'value': res_name,      'type': 'INT',     'domain': 'POINT'},
        {'name': 'atomic_number',   'value': atomic_number, 'type': 'INT',     'domain': 'POINT'},
        {'name': 'b_factor',        'value': mol.b_factor,  'type': 'FLOAT',   'domain': 'POINT'},
        {'name': 'vdw_radii',       'value': vdw_radii,     'type': 'FLOAT',   'domain': 'POINT'},
        {'name': 'chain_id',        'value': chain_id,      'type': 'INT',     'domain': 'POINT'},
        {'name': 'is_backbone',     'value': is_backbone,   'type': 'BOOLEAN', 'domain': 'POINT'},
        {'name': 'is_alpha_carbon', 'value': is_alpha,      'type': 'BOOLEAN', 'domain': 'POINT'},
        {'name': 'is_hetero',       'value': mol.hetero,    'type': 'BOOLEAN', 'domain': 'POINT'},
        {'name': 'is_solvent',      'value': is_solvent,    'type': 'BOOLEAN', 'domain': 'POINT'},
        {'name': 'is_nucleic',      'value': is_nucleic,    'type': 'BOOLEAN', 'domain': 'POINT'},
        {'name': 'is_peptide',      'value': is_peptide,    'type': 'BOOLEAN', 'domain': 'POINT'},
        {'name': 'is_carb',         'value': is_carb,       'type': 'BOOLEAN', 'domain': 'POINT'}
    )
    
    # assign the attributes to the object
    for att in attributes:
        add_attribute(mol_object, att['name'], att['value'], att['type'], att['domain'])
    
    # add custom properties for the object, such as number of chains, biological assemblies etc
    mol_object['chain_id_unique'] = list(np.unique(mol.chain_id))
    
    return mol_object
