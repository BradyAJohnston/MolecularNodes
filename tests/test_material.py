import bpy
from numpy.testing import assert_allclose

import molecularnodes as mn


def test_setting_material():
    mol = mn.Molecule.fetch("4ozs").add_style("cartoon")
    s = mol.styles[0]
    s.material = "MN Default"
    assert isinstance(s.material, mn.material.MaterialTreeInterface)
    assert isinstance(s.material.material, bpy.types.Material)
    s.material = mn.material.AmbientOcclusion()
    assert isinstance(s.material, mn.material.MaterialTreeInterface)
    assert isinstance(s.material.material, bpy.types.Material)


def test_generic_material():
    mol = mn.Molecule.fetch("4ozs").add_style("cartoon")
    s = mol.styles[0]
    s.material = bpy.data.materials["Material"]
    assert isinstance(s.material, mn.material.MaterialTreeInterface)
    assert isinstance(s.material.material, bpy.types.Material)

    assert_allclose(s.material.principled_bsdf_base_color, (0.8, 0.8, 0.8, 1.0))


def test_ambient_occlusion():
    mol = mn.Molecule.fetch("4ozs").add_style("cartoon")
    s = mol.styles[0]
    s.material = "MN Ambient Occlusion"
    assert isinstance(s.material, mn.material.MaterialTreeInterface)
    assert isinstance(s.material.material, bpy.types.Material)

    assert_allclose(s.material.ambient_occlusion_distance, 1.0)
    s.material.ambient_occlusion_distance = 0.1
    assert_allclose(s.material.ambient_occlusion_distance, 0.1)


def test_material_attribute_access():
    for material in [
        mn.material.AmbientOcclusion,
        mn.material.Default,
        mn.material.FlatOutline,
        mn.material.Squishy,
        mn.nodes.material.dynamic_material_interface(bpy.data.materials["Material"]),
    ]:
        for attr in material.__dict__:
            if not attr.startswith("_"):
                assert hasattr(material, attr)
                assert getattr(material, attr) is not None
