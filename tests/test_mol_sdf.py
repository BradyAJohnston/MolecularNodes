import molecularnodes as mn
import pytest
import bpy

from .constants import test_data_directory, attributes
from .utils import sample_attribute_to_string

mn.unregister()
mn.register()


def evaluate(object):
    object.update_tag()
    dg = bpy.context.evaluated_depsgraph_get()
    return object.evaluated_get(dg)


formats = ['mol', 'sdf']


@pytest.mark.parametrize("format", formats)
def test_open(snapshot, format):
    mol, file = mn.io.local.open_structure_local_mol(
        test_data_directory / f'caffeine.{format}')

    assert mol
    assert file


@pytest.mark.parametrize("format", formats)
@pytest.mark.parametrize("style", ['ball_and_stick', 'spheres', 'surface'])
def test_load(snapshot, format, style):
    mol = mn.io.local.load(test_data_directory /
                           f'caffeine.{format}', style=style)
    if style == 'spheres':
        mn.blender.nodes.get_style_node(
            mol).inputs['EEVEE'].default_value = True
    mn.blender.nodes.realize_instances(mol)

    for attribute in attributes:
        snapshot.assert_match(
            sample_attribute_to_string(evaluate(mol), attribute),
            f"{attribute}.txt"
        )
