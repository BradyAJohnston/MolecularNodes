import bpy
import pytest
from numpy.testing import assert_allclose
import molecularnodes as mn


def test_setting_material():
    mol = mn.Molecule.fetch("4ozs").add_style("cartoon")
    s = mol.styles[0]
    s.material = "MN Default"
    assert isinstance(s.material, mn.material.TreeInterface)
    assert isinstance(s.material.material, bpy.types.Material)
    s.material = mn.material.AmbientOcclusion()
    assert isinstance(s.material, mn.material.TreeInterface)
    assert isinstance(s.material.material, bpy.types.Material)
    with pytest.raises(AttributeError):
        s.material.non_existent_property
    with pytest.raises(AttributeError):
        s.material.non_existent_property = 1.0
    with pytest.raises(mn.nodes.interface.SocketLinkedError):
        s.material.mix_a


def test_setting_material_in_apply_builtins():
    materials = [
        mn.material.AmbientOcclusion(),
        mn.material.Default(),
        mn.material.FlatOutline(),
        mn.material.Squishy(),
        mn.material.TransparentOutline(),
    ]
    for material in materials:
        mol = mn.Molecule.fetch("4ozs").add_style("cartoon", material=material)
        assert mol is not None, f"Failed to apply {material.__class__.__name__}"


def test_setting_material_in_apply_custom():
    def create_new_material(name="New_Material", color=(1, 0, 0, 1)):
        if name in bpy.data.materials:
            material = bpy.data.materials[name]
        else:
            material = bpy.data.materials.new(name=name)
        material.use_nodes = True
        material.node_tree.nodes.clear()
        principled_node = material.node_tree.nodes.new("ShaderNodeBsdfPrincipled")
        principled_node.location = (0, 0)
        principled_node.inputs["Base Color"].default_value = color
        principled_node.inputs["Metallic"].default_value = 0.5
        principled_node.inputs["Roughness"].default_value = 0.2
        output_node = material.node_tree.nodes.new("ShaderNodeOutputMaterial")
        output_node.location = (300, 0)
        material.node_tree.links.new(
            principled_node.outputs["BSDF"], output_node.inputs["Surface"]
        )
        return material

    new_mat = create_new_material(name="My_Red_Material", color=(1, 0, 0, 1))
    mol = mn.Molecule.fetch("4ozs").add_style("cartoon", material=new_mat)
    assert mol is not None


def test_generic_material():
    mol = mn.Molecule.fetch("4ozs").add_style("cartoon")
    s = mol.styles[0]
    s.material = bpy.data.materials["Material"]
    assert s.material.material.name == "Material"  # type: ignore
    with pytest.raises(KeyError):
        s.material = "Non-existent material"
    with pytest.raises(TypeError):
        s.material = 1
    assert isinstance(s.material, mn.material.TreeInterface)
    assert isinstance(s.material.material, bpy.types.Material)

    assert_allclose(s.material.principled_bsdf_base_color, (0.8, 0.8, 0.8, 1.0))
    s.material = "MN Ambient Occlusion"
    assert s.material.material.name == "MN Ambient Occlusion"  # type: ignore


def test_ambient_occlusion():
    mol = mn.Molecule.fetch("4ozs").add_style("cartoon")
    s = mol.styles[0]
    assert isinstance(s.material, mn.material.TreeInterface)
    assert isinstance(s.material.material, bpy.types.Material)
    s.material = "MN Ambient Occlusion"
    assert s.material.material.name == "MN Ambient Occlusion"

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
