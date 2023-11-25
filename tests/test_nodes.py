import bpy
import pytest
import molecularnodes as mn
from molecularnodes.blender import nodes
def test_node_name_format():
    assert mn.blender.nodes.format_node_name("MN_style_cartoon") == "Style Cartoon"

def test_get_nodes():
    mol = mn.io.pdb.load('4ozs', style='cartoon')
    
    assert nodes.get_nodes_last_output(mol.modifiers['MolecularNodes'].node_group)[0].name == "MN_style_cartoon"
    assert nodes.get_style_node(mol).name == "MN_style_cartoon"
    
    mol2 = mn.io.pdb.load('1cd3', style='cartoon', build_assembly=True)
    
    assert nodes.get_nodes_last_output(mol2.modifiers['MolecularNodes'].node_group)[0].name == "MN_assembly_1cd3"
    assert nodes.get_style_node(mol2).name == "MN_style_cartoon"
