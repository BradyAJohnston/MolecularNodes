import MolecularNodes as mn
import bpy
import pytest

@pytest.fixture
def bpy_module(cache):
    return cache.get("bpy_module", None)

def test_node_surface(bpy_module):
    obj = mn.load.molecule_rcsb('6n2y')
    name = 'MOL_style_surface_split_6n2y'
    split_surface_node = mn.nodes.create_custom_surface(name, len(obj['chain_id_unique']))
    assert split_surface_node.name == name