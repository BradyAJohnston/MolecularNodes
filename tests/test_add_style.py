import bpy
import pytest
import molecularnodes as mn
from molecularnodes.nodes import node_management


def test_add_style_with_selection():
    mol = mn.Molecule.fetch("4ozs").add_style("cartoon")
    mol.store_named_attribute(mol.named_attribute("is_side_chain"), "show_side_chains")
    mol.add_style("ball_and_stick", selection="show_side_chains")

    node_style = [
        node
        for node in mol.modifier_node_tree
        if (
            isinstance(node, bpy.types.NodeGroup)
            and node.node_tree.name == "Style Ball and Stick"
        )
    ][0]

    assert (
        node_style.inputs["Selection"].links[0].from_node.inputs["Name"].default_value
        == "show_side_chains"
    )

    with pytest.warns(UserWarning):
        mol.add_style("cartoon", selection="non_existing_selection")
