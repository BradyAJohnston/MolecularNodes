import pytest
import molecularnodes as mn
from molecularnodes.nodes import node_management


def test_add_style_with_selection():
    mol = mn.Molecule.fetch("4ozs").add_style("cartoon")
    mol.store_named_attribute(mol.named_attribute("is_side_chain"), "show_side_chains")
    mol.add_style("ball_and_stick", selection="show_side_chains")

    node_style = mol.modifier_node_tree.nodes["Style Ball and Stick"]
    assert (
        node_style.inputs["Selection"].links[0].from_node.inputs["Name"].default_value
        == "show_side_chains"
    )

    with pytest.warns(UserWarning):
        mol.add_style("cartoon", selection="non_existing_selection")


def test_change_style_values():
    mol = mn.Molecule.fetch("4ozs").add_style("cartoon")
    pre = mol.named_attribute("position", evaluate=True)
    style_node = node_management.get_final_style_nodes(mol.modifier_node_tree)[0]
    style_node.inputs["Quality"].default_value = 5
    post = mol.named_attribute("position", evaluate=True)

    assert len(pre) < len(post)
