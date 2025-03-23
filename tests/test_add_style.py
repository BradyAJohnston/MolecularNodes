import molecularnodes as mn
from numpy.testing import assert_allclose
import numpy as np
import bpy
import pytest


def test_style_interface():
    mol = mn.Molecule.fetch("4ozs").add_style("cartoon")
    assert len(mol.tree.nodes) == 7
    mn.blender.styles.add_style_branch(mol.tree, "cartoon")
    assert len(mol.tree.nodes) == 8
    mn.blender.styles.input_named_attribute(
        mol.tree.nodes["Style Cartoon"].inputs["Selection"], "is_backbone", "BOOLEAN"
    )
    assert len(mol.tree.nodes) == 9
    mn.blender.styles.add_style_branch(mol.tree, "surface", selection="is_backbone")
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
    mn.blender.styles.add_style_branch(mol.tree, "spheres")

    # testing the current interface for node trees via scripting. We can
    # expose their values through helper classes
    w = mn.blender.styles.StyleWrangler(mol.tree)
    style = w.styles[3]
    assert_allclose(
        style.cartoon_width,
        bpy.data.node_groups["Style Cartoon"]
        .interface.items_tree["Width"]
        .default_value,
    )
    style.cartoon_width = 1.0
    assert_allclose(style.cartoon_width, 1.0)


def test_add_color_node():
    mol = mn.Molecule.fetch("4ozs").add_style("spheres")
    assert len(mol.tree.nodes) == 7
    mn.blender.styles.add_style_branch(mol.tree, "spheres")
    assert len(mol.tree.nodes) == 8
    # if we are adding a style with a Set Color node, we check that 3 extra nodes
    # have been added rather than just 1 (style, color & named attribute), then we check
    # that the Set Color nodes has an input for the "Color" socket that is a named attribute
    # node, checking that the name is the one that we set
    mn.blender.styles.add_style_branch(mol.tree, "cartoon", color="is_peptide")
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

    with pytest.raises(UserWarning):
        mol.add_style("cartoon", selection="non_existing_selection")
