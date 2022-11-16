import bpy
import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdb as pdb
import biotite.structure.io.pdbx as pdbx
import biotite.structure.io.mmtf as mmtf
import biotite.database.rcsb as rcsb

def open_structure_rcsb(pdb_code, include_bonds = True):
    file = mmtf.MMTFFile.read(rcsb.fetch(pdb_code, "mmtf"))
    mol = mmtf.get_structure(file, model = 1, extra_fields = ["b_factor", "charge"], include_bonds = include_bonds)
    return mol


def open_structure_local_pdb(file_path, include_bonds = True):
    file = pdb.PDBFile.read(file_path)
    mol = pdb.get_structure(file, model = 1, extra_fields = ['b_factor', 'charge'], include_bonds = include_bonds)
    return mol

def open_structure_local_pdbx(file_path, include_bonds = True):
    file = pdbx.PDBxFile.read(file_path)
    mol  = pdbx.get_structure(file, model = 1, extra_fields = ['b_factor', 'charge'])
    # for some reason .pdbx doesn't have an option for bonds on import, so manually create
    # them here if requested
    if include_bonds:
        mol.bonds = struc.bonds.connect_via_residue_names(mol, inter_residue = True)
    return mol

def mn_collection():
    coll = bpy.data.collections.get('MolecularNodes')
    if not coll:
        coll = bpy.data.collections.new('MolecularNodes')
        bpy.context.scene.collection.children.link(coll)
    return coll

def create_model(name, collection, locations, bonds=[]):
    """
    Creates a mesh with the given name in the given collection, from the supplied
    values for the locationso of vertices, and if supplied, bonds and faces.
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
    attribute = object.data.attributes.new(name, type, domain)
    attribute.data.foreach_set('value', data)

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
    centroid = struc.centroid(mol)* world_scale

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
    chain_id = np.searchsorted(np.unique(mol.chain_id), mol.chain_id) + 1 # add 1 to start from chain 1
    vdw_radii =  np.fromiter(map(struc.info.vdw_radius_single, mol.element), dtype=np.float) * world_scale
    is_alpha = np.fromiter(map(lambda x: x == "CA", mol.atom_name), dtype = np.bool)
    is_solvent = struc.filter_solvent(mol)

    if include_bonds:
        bond_type = bonds[:, 2].copy(order = 'C') # the .copy(order = 'C') is to fix a weird ordering issue with the resulting array
        add_attribute(mol_object, 'bond_type', bond_type,"INT", "EDGE")

    # assign the attributes to the object
    add_attribute(mol_object, 'res_id', mol.res_id, "INT")
    add_attribute(mol_object, 'res_name', res_name, "INT")
    add_attribute(mol_object, "atomic_number", atomic_number, "INT", "POINT")
    add_attribute(mol_object, "b_factor", mol.b_factor, "FLOAT", "POINT")
    add_attribute(mol_object, "is_backbone", struc.filter_backbone(mol),"BOOLEAN", "POINT")
    add_attribute(mol_object, "is_alpha_carbon", is_alpha, "BOOLEAN", "POINT")
    add_attribute(mol_object, 'is_hetero', mol.hetero, "BOOLEAN")
    add_attribute(mol_object, "vdw_radii", vdw_radii, "FLOAT", "POINT")
    add_attribute(mol_object, "chain_id", chain_id, "INT", "POINT")
    add_attribute(mol_object, 'is_solvent', is_solvent, "BOOLEAN", "POINT")