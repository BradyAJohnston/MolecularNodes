import bpy
import pytest
import molecularnodes as mn
from molecularnodes.nodes import material
from .constants import data_dir

PRESETS = [
    material.Default,
    material.AmbientOcclusion,
    material.Flat,
    material.Squishy,
    material.TransparentOutline,
]


@pytest.mark.parametrize("preset", PRESETS)
def test_preset_creates_material(preset):
    mat = preset()
    assert isinstance(mat.material, bpy.types.Material)
    assert mat.material.name.startswith(preset.name)
    assert len(mat.material.node_tree.nodes) > 0


def test_presets_are_independent():
    a = material.AmbientOcclusion()
    b = material.AmbientOcclusion()
    assert a.material != b.material
    a.distance = 5.0
    assert b.distance == pytest.approx(1.0)


def test_preset_parameters_map_to_sockets():
    ao = material.AmbientOcclusion(distance=0.15, exponent=3.0)
    assert ao.distance == pytest.approx(0.15)
    assert ao.exponent == pytest.approx(3.0)
    assert ao.ao.i.distance.default_value == pytest.approx(0.15)

    ao.distance = 2.0
    assert ao.ao.i.distance.default_value == pytest.approx(2.0)


def test_flat_outline_toggle():
    flat = material.Flat(outline=False, threshold=0.5, thickness=0.3)
    assert flat.outline is False
    assert flat.node.i.outline.default_value == "None"
    flat.outline = True
    assert flat.node.i.outline.default_value == "Outline"
    assert flat.threshold == pytest.approx(0.5)
    assert flat.thickness == pytest.approx(0.3)


def test_transparent_outline():
    mat = material.TransparentOutline(alpha=0.5, outline=False)
    assert mat.alpha == pytest.approx(0.5)
    assert mat.outline is False
    assert mat.node.i.menu.default_value == "Transparent"
    mat.outline = True
    assert mat.node.i.menu.default_value == "Outline"
    assert mat.material.surface_render_method == "BLENDED"


def test_default_material_parameters():
    mat = material.Default(roughness=0.5, ao_distance=1.0)
    assert mat.roughness == pytest.approx(0.5)
    assert mat.ao_distance == pytest.approx(1.0)
    assert mat.bsdf.i.roughness.default_value == pytest.approx(0.5)


def test_squishy_material_parameters():
    mat = material.Squishy(subsurface_scale=0.5)
    assert mat.subsurface_scale == pytest.approx(0.5)
    assert mat.bsdf.i.subsurface_weight.default_value == pytest.approx(1.0)


def test_add_style_with_preset():
    mol = mn.Molecule.fetch("1BNA", cache=data_dir)
    mat = material.AmbientOcclusion(distance=0.5)
    mol.add_style("spheres", material=mat)
    socket = mol.tree.tree.nodes["Style Spheres"].inputs["Material"]
    assert socket.default_value == mat.material

    mat.distance = 3.0
    assert mat.ao.i.distance.default_value == pytest.approx(3.0)


def test_add_style_with_material_string():
    mol = mn.Molecule.fetch("1BNA", cache=data_dir)
    mol.add_style("cartoon", material="MN Flat")
    socket = mol.tree.tree.nodes["Style Cartoon"].inputs["Material"]
    assert socket.default_value.name.startswith("MN Flat")


def test_append_material():
    mat = material.append_material("MN Default")
    assert isinstance(mat, bpy.types.Material)
