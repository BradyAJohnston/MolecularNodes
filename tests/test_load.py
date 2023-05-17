import MolecularNodes as mn
import pytest
import os
import bpy

# ensure we can successfully install all of the required pacakges 
# def test_install_packages():
    # mn.pkg.install_all_packages()
    # assert mn.pkg.is_current('biotite') == True

def get_verts(obj, float_decimals = 4, n_verts = 100, apply_modifiers = True):
    """
    Get the first n_verts number of verts from an object.
    """
    if apply_modifiers:
        bpy.context.view_layer.objects.active = obj
        for modifier in obj.modifiers:
            bpy.ops.object.modifier_apply(modifier = modifier.name)
    verts = ""
    for i, v in enumerate(obj.data.vertices):
        if i < n_verts:
            vert = [v.co.x, v.co.y, v.co.z]
            vert = list(map(lambda x: round(x, float_decimals), vert))
            verts += "{},{},{}\n".format(vert[0], vert[1], vert[2])
    return verts

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

def test_local_pdb(snapshot):
    files = [f"tests/data/1l58.{ext}" for ext in ['cif', 'pdb']]
    obj1, obj2 = map(mn.load.molecule_local, files)
    obj3 = mn.load.molecule_rcsb('1l58')
    verts_1, verts_2, verts_3 = map(lambda x: get_verts(x, apply_modifiers = False), [obj1, obj2, obj3])
    assert verts_1 == verts_2
    assert verts_1 == verts_3
    snapshot.assert_match(verts_1, '1L58_verts.txt')