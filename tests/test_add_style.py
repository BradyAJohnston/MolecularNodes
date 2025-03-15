import molecularnodes as mn
from numpy.testing import assert_allclose
import bpy


def test_style_interface():
    mol = mn.Molecule.fetch("4ozs").add_style("cartoon")
    assert len(mol.tree.nodes) == 6
    mn.blender.styles.add_style_branch(mol.tree, "cartoon")
    assert len(mol.tree.nodes) == 8
    mn.blender.styles.add_style_branch(mol.tree, "surface")
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
