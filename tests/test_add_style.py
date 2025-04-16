import bpy
import numpy as np
import pytest
from numpy.testing import assert_allclose
import molecularnodes as mn
from molecularnodes.nodes import geometry
from molecularnodes.nodes.geometry import add_style_branch


def test_style_interface():
    mol = mn.Molecule.fetch("4ozs").add_style("cartoon")
    assert len(mol.tree.nodes) == 7
    add_style_branch(mol.tree, "cartoon")
    assert len(mol.tree.nodes) == 8
    geometry.input_named_attribute(
        mol.tree.nodes["Style Cartoon"].inputs["Selection"], "is_backbone", "BOOLEAN"
    )
    assert len(mol.tree.nodes) == 9
    add_style_branch(mol.tree, "surface", selection="is_backbone")
    assert len(mol.tree.nodes) == 11
    print(f"{list(mol.tree.nodes)}")

    assert (
        mol.tree.nodes["Style Surface"]
        .inputs["Selection"]
        .links[0]
        .from_node.inputs["Name"]
        .default_value  # type: ignore
        == "is_backbone"
    )
    add_style_branch(mol.tree, "spheres")

    style = mol.styles[0]
    assert_allclose(
        style.peptide_width,
        bpy.data.node_groups["Style Cartoon"]
        .interface.items_tree["Peptide Width"]
        .default_value,
    )
    style.peptide_width = 1.0
    assert_allclose(style.peptide_width, 1.0)

    assert len(mol.tree.nodes) == 12
    mol.add_style("cartoon", color="is_peptide")
    assert len(mol.tree.nodes) == 15
    mol.styles[-1].remove()
    assert len(mol.tree.nodes) == 12


def test_add_color_node():
    mol = mn.Molecule.fetch("4ozs").add_style("spheres")
    assert len(mol.tree.nodes) == 7
    add_style_branch(mol.tree, "spheres")
    assert len(mol.tree.nodes) == 8
    # if we are adding a style with a Set Color node, we check that 3 extra nodes
    # have been added rather than just 1 (style, color & named attribute), then we check
    # that the Set Color nodes has an input for the "Color" socket that is a named attribute
    # node, checking that the name is the one that we set
    add_style_branch(mol.tree, "cartoon", color="is_peptide")
    assert len(mol.tree.nodes) == 11
    node_sc = mol.tree.nodes["Style Cartoon"].inputs[0].links[0].from_node  # type: ignore
    assert node_sc.inputs["Color"].is_linked
    node_na = node_sc.inputs["Color"].links[0].from_socket.node
    assert node_na.inputs["Name"].default_value == "is_peptide"
    assert node_na.data_type == "FLOAT_COLOR"


def test_add_style_with_selection():
    mol = mn.Molecule.fetch("4ozs").add_style("cartoon")
    mol.select.res_id(range(50)).is_side_chain().store_selection("show_side_chains")
    mol.add_style("ball+stick", selection="show_side_chains")
    mol.add_style("cartoon")

    sel = (
        mn.entities.MoleculeSelector()
        .res_id(range(50, 500))
        .is_side_chain()
        .res_name(["ARG", "LYS", "VAL"])
    )

    mol.add_style(
        "ball+stick",
        selection=sel,
    )

    assert "sel_0" in mol.list_attributes()
    assert np.allclose(mol.named_attribute("sel_0"), sel.evaluate_on_array(mol.array))

    with pytest.warns(UserWarning):
        mol.add_style("cartoon", selection="non_existing_selection")


def test_change_style_values():
    mol = mn.Molecule.fetch("4ozs").add_style("cartoon")
    pre = mol.named_attribute("position", evaluate=True)
    mol.styles[0].quality = 5
    post = mol.named_attribute("position", evaluate=True)

    assert len(pre) < len(post)
    with pytest.raises(ValueError):
        mol.styles[0].quality = 1.0
