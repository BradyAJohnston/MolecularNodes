import molecularnodes as mn
from numpy.testing import assert_allclose
import bpy


def test_style_interface():
    mol = mn.Molecule.fetch("4ozs").add_style("cartoon")
    assert len(mol.tree.nodes) == 6
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
    assert len(mol.tree.nodes) == 6

    mn.blender.styles.add_style_branch(mol.tree, "cartoon", color="position")
    assert len(mol.tree.nodes) == 10
    node_sc = mol.tree.nodes["Style Cartoon"].inputs[0].links[0].from_node
    assert node_sc.inputs["Color"].is_linked
    print(f"{list(node_sc.inputs)}")
    assert (
        node_sc.inputs["Color"].links[0].from_socket.node.inputs["Name"].default_value
        == "Position"
    )
