import bpy
from numpy.testing import assert_allclose

import molecularnodes as mn


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

    assert mol.material.power == 1.5  # type: ignore
    mol.material.power = 0.1
    assert_allclose(mol.material.power, 0.1)
    assert mol.material.samples == 16
    mol.material.samples = 3
    assert mol.material.samples == 3


def test_material_attribute_access():
    for material in [
        mn.material.AmbientOcclusion,
        mn.material.Default,
        mn.material.FlatOutline,
        mn.material.Squishy,
        mn.material.generic_material_interface(bpy.data.materials["Material"]),
    ]:
        for attr in material.__dict__:
            if not attr.startswith("_"):
                assert hasattr(material, attr)
                assert getattr(material, attr) is not None
