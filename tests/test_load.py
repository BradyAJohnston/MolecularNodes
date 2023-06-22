import bpy
import os
import pytest
import MolecularNodes as mn

# ensure we can successfully install all of the required pacakges 
# def test_install_packages():
    # mn.pkg.install_all_packages()
    # assert mn.pkg.is_current('biotite') == True

def apply_mods(obj):
    """
    Applies the modifiers on the modifier stack
    
    This will realise the computations inside of any Geometry Nodes modifiers, ensuring
    that the result of the node trees can be compared by looking at the resulting 
    vertices of the object.
    """
    bpy.context.view_layer.objects.active = obj
    for modifier in obj.modifiers:
        bpy.ops.object.modifier_apply(modifier = modifier.name)

def get_verts(obj, float_decimals=4, n_verts=100, apply_modifiers=True, seed=42):
    """
    Randomly samples a specified number of vertices from an object.

    Parameters
    ----------
    obj : object
        Object from which to sample vertices.
    float_decimals : int, optional
        Number of decimal places to round the vertex coordinates, defaults to 4.
    n_verts : int, optional
        Number of vertices to sample, defaults to 100.
    apply_modifiers : bool, optional
        Whether to apply all modifiers on the object before sampling vertices, defaults to True.
    seed : int, optional
        Seed for the random number generator, defaults to 42.

    Returns
    -------
    str
        String representation of the randomly selected vertices.

    Notes
    -----
    This function randomly samples a specified number of vertices from the given object.
    By default, it applies all modifiers on the object before sampling vertices. The
    random seed can be set externally for reproducibility.

    If the number of vertices to sample (`n_verts`) exceeds the number of vertices
    available in the object, all available vertices will be sampled.

    The vertex coordinates are rounded to the specified number of decimal places
    (`float_decimals`) before being included in the output string.

    Examples
    --------
    >>> obj = mn.load.molecule_rcsb('6n2y', starting_style=2)
    >>> get_verts(obj, float_decimals=3, n_verts=50, apply_modifiers=True, seed=42)
    '1.234,2.345,3.456\n4.567,5.678,6.789\n...'
    """

    import random

    random.seed(seed)

    if apply_modifiers:
        apply_mods(obj)

    vert_list = [(v.co.x, v.co.y, v.co.z) for v in obj.data.vertices]

    if n_verts > len(vert_list):
        n_verts = len(vert_list)

    random_verts = random.sample(vert_list, n_verts)

    verts_string = ""
    for i, vert in enumerate(random_verts):
        if i < n_verts:
            rounded = [round(x, float_decimals) for x in vert]
            verts_string += "{},{},{}\n".format(rounded[0], rounded[1], rounded[2])

    return verts_string


def test_open_rcsb(snapshot):
    mn.load.open_structure_rcsb('4ozs')
    assert True == True

def test_rcsb_4ozs(snapshot):
    obj = mn.load.molecule_rcsb('4ozs')
    verts = get_verts(obj, apply_modifiers = False)
    snapshot.assert_match(verts, '4ozs_verts.txt')

def test_rcsb_6n2y_cartoon(snapshot):
    obj = mn.load.molecule_rcsb('6n2y', starting_style=2)
    verts = get_verts(obj)
    snapshot.assert_match(verts, '6n2y_cartoon_verts.txt')

def test_rcsb_6n2y_ribbon(snapshot):
    obj = mn.load.molecule_rcsb('6n2y', starting_style=3)
    verts = get_verts(obj)
    snapshot.assert_match(verts, '6n2y_ribbon_verts.txt')

def test_rcsb_6n2y_surface_split(snapshot):
    obj = mn.load.molecule_rcsb('6n2y', starting_style=1, setup_nodes = True)
    node_surface = mn.nodes.create_custom_surface(
        name = 'MOL_style_surface_6n2y_split', 
        n_chains = len(obj['chain_id_unique'])
        )
    node_group = obj.modifiers['MolecularNodes'].node_group
    node_group.nodes['Group.001'].node_tree = node_surface
    
    for link in node_group.links:
        if link.to_node.name == "Group.001":
            node_group.links.remove(link)
    new_link = node_group.links.new
    new_link(
        node_group.nodes['Group'].outputs[0], 
        node_group.nodes['Group.001'].inputs[0]
    )
    new_link(
        node_group.nodes['Group.001'].outputs[0], 
        node_group.nodes['Group Output'].inputs[0]
    )
    
    verts = get_verts(obj, n_verts=1000, float_decimals=2)
    snapshot.assert_match(verts, '6n2y_surface_verts.txt')

def test_local_pdb(snapshot):
    files = [f"tests/data/1l58.{ext}" for ext in ['cif', 'pdb']]
    obj1, obj2 = map(mn.load.molecule_local, files)
    obj3 = mn.load.molecule_rcsb('1l58')
    verts_1, verts_2, verts_3 = map(lambda x: get_verts(x, apply_modifiers = False), [obj1, obj2, obj3])
    assert verts_1 == verts_2
    assert verts_1 == verts_3
    snapshot.assert_match(verts_1, '1L58_verts.txt')

def test_esmfold(snapshot):
    sequence = "HHHHHH"
    obj = mn.esmfold.molecule_esmfold(sequence)
    verts = get_verts(obj, apply_modifiers = False)
    snapshot.assert_match(verts, 'esmfold_6xHis_verts.txt')

def test_starfile_positions(snapshot):
    file = "tests/data/cistem.star"
    obj = mn.star.load_star_file(file)
    verts = get_verts(obj, n_verts = 500, apply_modifiers = False)
    snapshot.assert_match(verts, 'starfile_verts.txt')

def test_md_load_gro_xtc(snapshot):
    top = "tests/data/md_ppr/box.gro"
    traj = "tests/data/md_ppr/first_5_frames.xtc"
    
    obj, coll = mn.md.load_trajectory(top, traj)
    verts = get_verts(obj, apply_modifiers = False)
    
    snapshot.assert_match(verts, 'md_gro_xtc_verts.txt')

def test_rcsb_nmr(snapshot):
    CODE = "2M6Q"
    obj = mn.load.molecule_rcsb(CODE)
    coll_frames = bpy.data.collections.get(CODE + "_frames")
    assert len(coll_frames.objects) == 10
    
    verts = get_verts(obj, apply_modifiers = False)
    snapshot.assert_match(verts, 'rcsb_nmr_2M6Q.txt')

def test_load_small_mol(snapshot):
    file = "tests/data/ASN.cif"
    obj = mn.load.molecule_local(file)
    verts = get_verts(obj, apply_modifiers = False)
    snapshot.assert_match(verts, 'asn_atoms.txt')
    
    bond_types = mn.obj.get_attribute(obj, 'bond_type')
    edges = ''.join([str(bond_type) for bond_type in bond_types])
    snapshot.assert_match(edges, 'asn_edges.txt')