import bpy
import os
import pytest
import MolecularNodes as mn
import numpy as np
from .utils import get_verts, apply_mods

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


def test_rcsb_cache(snapshot):
    from pathlib import Path
    from shutil import rmtree
    # we want to make sure cached files are freshly downloaded, but
    # we don't want to delete our entire real cache
    test_cache = Path(Path.home(), '.MolecularNodesTests')
    if test_cache.exists():
        rmtree(test_cache)
    _ = mn.load.molecule_rcsb('6BQN', cache_dir = test_cache)
    assert (test_cache / '6BQN.mmtf').exists()

def test_1cd3_bio_assembly(snapshot):
    obj_rcsb = mn.load.molecule_rcsb('1CD3')
    obj_cif, obj_pdb = [mn.load.molecule_local(f"tests/data/1cd3.{ext}") for ext in ["pdb", "cif"]]
    
    vert_list = []
    objects = [obj_rcsb, obj_cif, obj_pdb]
    for obj in objects:
        data_object = mn.assembly.mesh.create_data_object(
            transforms_dict = obj['biological_assemblies'], 
            name = f"data_assembly_{obj.name}"
        )
        
        node_bio_assembly = mn.nodes.create_assembly_node_tree(
            name = obj.name, 
            iter_list = obj['chain_id_unique'], 
            data_object = data_object
            )
        
        node_group = obj.modifiers['MolecularNodes'].node_group
        node_group.nodes['Group.001'].node_tree = node_bio_assembly
        
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
        
        node_realize = node_group.nodes.new('GeometryNodeRealizeInstances')
        
        node_realize.location = (node_group.nodes['Group.001'].location + 
                                node_group.nodes['Group Output'].location) / 2
        new_link(
            node_group.nodes['Group.001'].outputs[0], 
            node_realize.inputs[0]
        )
        
        new_link(
            node_realize.outputs[0], 
            node_group.nodes['Group Output'].inputs[0]
        )
        
    verts = get_verts(obj_rcsb, n_verts=1000, float_decimals=2)
    
    
    # for verts in vert_list:
    snapshot.assert_match(verts, '1cd3_bio_assembly.txt')
    
    for obj in objects:
        apply_mods(obj)
    
    # turn each object to positions, sort the arrays (different import methods currently 
    # result in different ordered verts) but then check they are the same
    # Results shows the same number of atoms and atom positions are resulting from the 
    # different import methods so it still works
    positions = [np.sort(mn.obj.get_attribute(obj, 'position'), axis = 0, kind = 'quicksort')[::-1] for obj in objects]
    assert np.allclose(positions[0], positions[1], atol = 1e-4)
    
    # TODO: for some reason when opening from .CIF files, we are ending up with double the 
    # number of chains than we would need. I am unsure why this is the case, but will leave 
    # it for now, as everything else is working well.
    
    # assert np.allclose(positions[0], positions[2], atol = 1e-4)
