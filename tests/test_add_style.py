import molecularnodes as mn
from numpy.testing import assert_allclose


def test_add_style():
    mol = mn.Molecule.fetch("4ozs").add_style("cartoon")

    assert len(mol.tree.nodes) == 6

    mn.blender.styles.add_style_branch(mol.tree, "surface")

    assert len(mol.tree.nodes) == 8


def test_style_interface():
    mol = mn.Molecule.fetch("4ozs").add_style("cartoon")

    mn.blender.styles.add_style_branch(mol.tree, "cartoon")
    mn.blender.styles.add_style_branch(mol.tree, "surface")
    mn.blender.styles.add_style_branch(mol.tree, "spheres")

    w = mn.blender.styles.StyleWrangler(mol.tree)
    style = w.styles[2]
    assert_allclose(style.cartoon_width, 2.2)
