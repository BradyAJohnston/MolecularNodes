import molecularnodes as mn
import bpy
from numpy.testing import assert_array_equal, assert_allclose


def test_generic_material():
    mol = mn.Molecule.fetch("4ozs").add_style("cartoon")
    mol.material = bpy.data.materials["Material"]
    assert isinstance(mol.material, mn.material.MaterialTreeInterface)
    assert isinstance(mol.material.material, bpy.types.Material)

    assert_allclose(mol.material.principled_bsdf_base_color, (0.8, 0.8, 0.8, 1.0))


def test_ambient_occlusion():
    mol = mn.Molecule.fetch("4ozs").add_style("cartoon")
    mol.material = "MN Ambient Occlusion"
    assert isinstance(mol.material, mn.material.MaterialTreeInterface)
    assert isinstance(mol.material.material, bpy.types.Material)

    assert mol.material.power == 1.5
    mol.material.power = 0.1
    assert_allclose(mol.material.power, 0.1)
