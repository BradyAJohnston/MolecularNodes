import MolecularNodes as mn
import bpy
import pytest

def get_verts(obj):
    verts = ""
    for v in obj.data.vertices:
        verts += "{},{},{}\n".format(v.co.x, v.co.y, v.co.z)
    return verts

def test_rcsb_4ozs(snapshot):
    obj = mn.load.molecule_rcsb('4ozs')
    snapshot.assert_match(get_verts(obj), '4ozs_verts.txt')

def test_rcsb_6n2y_cartoon(snapshot):
    obj = mn.load.molecule_rcsb('6n2y', starting_style=2)
    bpy.ops.object.modifier_apply(modifier="MolecularNodes")
    snapshot.assert_match(get_verts(obj), '6n2y_ribbon_verts.txt')

def test_rcsb_6n2y_ribbon(snapshot):
    obj = mn.load.molecule_rcsb('6n2y', starting_style=3)
    bpy.ops.object.modifier_apply(modifier="MolecularNodes")
    snapshot.assert_match(get_verts(obj), '6n2y_ribbon_verts.txt')
