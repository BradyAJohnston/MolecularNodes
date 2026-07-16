import pytest
import molecularnodes as mn
from molecularnodes.nodes import node_management
from molecularnodes.nodes.node_management import add_style_branch


def test_style_interface():
    mol = mn.Molecule.fetch("4ozs").add_style("cartoon")
    assert len(mol.modifier_node_tree.nodes) == 7
    add_style_branch(mol.modifier_node_tree, "cartoon")
    assert len(mol.modifier_node_tree.nodes) == 8
    node_management.input_named_attribute(
        mol.modifier_node_tree.nodes["Style Cartoon"].inputs["Selection"],
        "is_backbone",
        "BOOLEAN",
    )
    assert len(mol.modifier_node_tree.nodes) == 9
    add_style_branch(mol.modifier_node_tree, "surface", selection="is_backbone")
    assert len(mol.modifier_node_tree.nodes) == 11
    print(f"{list(mol.modifier_node_tree.nodes)}")

    assert (
        mol.modifier_node_tree.nodes["Style Surface"]
        .inputs["Selection"]
        .links[0]
        .from_node.inputs["Name"]
        .default_value
        == "is_backbone"
    )
    add_style_branch(mol.modifier_node_tree, "spheres")

    assert len(mol.modifier_node_tree.nodes) == 12
    mol.add_style("cartoon", color="is_peptide")
    assert len(mol.modifier_node_tree.nodes) == 15
    node_management.remove_style_node(
        node_management.get_final_style_nodes(mol.modifier_node_tree)[-1]
    )
    assert len(mol.modifier_node_tree.nodes) == 12


def test_add_color_node():
    mol = mn.Molecule.fetch("4ozs").add_style("spheres")
    assert len(mol.modifier_node_tree.nodes) == 7
    add_style_branch(mol.modifier_node_tree, "spheres")
    assert len(mol.modifier_node_tree.nodes) == 8
    # if we are adding a style with a Set Color node, we check that 3 extra nodes
    # have been added rather than just 1 (style, color & named attribute), then we check
    # that the Set Color nodes has an input for the "Color" socket that is a named attribute
    # node, checking that the name is the one that we set
    add_style_branch(mol.modifier_node_tree, "cartoon", color="is_peptide")
    assert len(mol.modifier_node_tree.nodes) == 11
    node_sc = mol.modifier_node_tree.nodes["Style Cartoon"].inputs[0].links[0].from_node
    assert node_sc.inputs["Color"].is_linked
    node_na = node_sc.inputs["Color"].links[0].from_socket.node
    assert node_na.inputs["Name"].default_value == "is_peptide"
    assert node_na.data_type == "FLOAT_COLOR"


def test_add_style_with_selection():
    mol = mn.Molecule.fetch("4ozs").add_style("cartoon")
    mol.store_named_attribute(mol.named_attribute("is_side_chain"), "show_side_chains")
    mol.add_style("ball+stick", selection="show_side_chains")

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
