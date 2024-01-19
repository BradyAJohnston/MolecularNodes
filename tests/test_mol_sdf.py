import molecularnodes as mn
import molecularnodes.blender as bl
import pytest
import bpy

from .constants import test_data_directory, attributes
from .utils import sample_attribute_to_string

mn.unregister()
mn.register()


formats = ['mol', 'sdf']


@pytest.mark.parametrize("format", formats)
def test_open(snapshot, format):
    molecule = mn.io.parse.SDF(test_data_directory / f'caffeine.{format}')

    assert molecule.array
    assert molecule.file


@pytest.mark.parametrize("format", formats)
@pytest.mark.parametrize("style", ['ball_and_stick', 'spheres', 'surface'])
def test_load(snapshot, format, style):
    object = mn.io.local.load(test_data_directory /
                              f'caffeine.{format}', style=style)
    if style == 'spheres':
        bl.nodes.get_style_node(object).inputs['EEVEE'].default_value = True
    mn.blender.nodes.realize_instances(object)

    for attribute in attributes:
        snapshot.assert_match(
            sample_attribute_to_string(bl.obj.evaluate(object), attribute),
            f"{attribute}.txt"
        )
