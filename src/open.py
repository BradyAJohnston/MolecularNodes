
def open_structure_rcsb(pdb_code):
    import biotite.structure as struc
    import biotite.database.rcsb as rcsb
    import biotite.structure.io.mmtf as mmtf

    file = mmtf.MMTFFile.read(rcsb.fetch(pdb_code, "mmtf"))
    mol = mmtf.get_structure(file, model = 1, extra_fields = ["b_factor", "charge"], include_bonds = True)
    return mol


def open_structure_local_pdb(file_path):
    import biotite.structure as struc
    import biotite.structure.io.pdb as pdb

    file = pdb.PDBFile.read(file_path)
    mol = pdb.get_structure(file, extra_fields = ['b_factor'], include_bonds = True)

    return mol


def import_protein_pdb(pdb, mol_name, center_molecule = False, del_solvent = False, include_bonds = True):
    import bpy
    import numpy as np
    import biotite.structure as struc
    import biotite.database.rcsb as rcsb
    import biotite.structure.io.mmtf as mmtf
    from . import packages
    from . import dict

    dict_elements = dict.dict_elements

    packages.install_packages()

    packages.verify()


    # remove the waters from the structure (most people don't want them anyway)
    if del_solvent:
        pdb = pdb[np.invert(struc.filter_solvent(pdb))]

    def mn_collection():
        
        coll = bpy.data.collections.get('MolecularNodes')
        if not coll:
            coll = bpy.data.collections.new('MolecularNodes')
            bpy.context.scene.collection.children.link(coll)
        return coll


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

    def add_attribute(object, name, data, type = "VECTOR", domain = "POINT"):
        attribute = object.data.attributes.new(name, type, domain)
        attribute.data.foreach_set('value', data)

    bonds = pdb.bonds.as_array()
    world_scale = 0.01
    locations = pdb.coord * world_scale
    centroid = struc.centroid(pdb)* world_scale

    # subtract the centroid from all of the positions to localise the molecule on the world origin
    if center_molecule:
        locations = locations - centroid

    if include_bonds:
        mod = create_model(
            name = mol_name, 
            collection = mn_collection(), 
            locations = locations, 
            bonds = bonds[:, [0,1]]
            )
        # set up relevant bond information attributes
        bond_type = bonds[:, 2].copy(order = 'C')
        add_attribute(mod, 'bond_type', bond_type,"INT", "EDGE")
    else:
        mod = create_model(
            name = mol_name, 
            collection = mn_collection(), 
            locations = locations
            )

    # compute some of the attributes
    is_solvent = struc.filter_solvent(pdb)
    chain_id = np.searchsorted(np.unique(pdb.chain_id), pdb.chain_id) + 1 # add 1 to start from chain 1
    atomic_number = np.fromiter(map(lambda x: dict_elements[x]["atomic_number"], np.char.title(pdb.element)), dtype = np.int)
    vdw_radii =  np.fromiter(map(struc.info.vdw_radius_single, pdb.element), dtype=np.float) * world_scale
    is_alpha = np.fromiter(map(lambda x: x == "CA", pdb.atom_name), dtype = np.bool)

    # assign the attributes to the object
    add_attribute(mod, 'res_id', pdb.res_id, "INT")
    add_attribute(mod, 'hetero', pdb.hetero, "BOOLEAN")
    add_attribute(mod, "b_factor", pdb.b_factor, "FLOAT", "POINT")
    add_attribute(mod, "is_alpha_carbon", is_alpha, "BOOLEAN", "POINT")
    add_attribute(mod, "atomic_number", atomic_number, "INT", "POINT")
    add_attribute(mod, "is_backbone", struc.filter_backbone(pdb),"BOOLEAN", "POINT")
    add_attribute(mod, "vdw_radii", vdw_radii,"FLOAT", "POINT")
    add_attribute(mod, "chain_id", chain_id, "INT", "POINT")
    add_attribute(mod, 'is_solvent', is_solvent, "BOOLEAN", "POINT")