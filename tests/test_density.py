import molecularnodes as mn
import bpy
import pytest
import numpy as np
import itertools
from .constants import data_dir
from .utils import (
    sample_attribute_to_string
)
try:
    import pyopenvdb
except ImportError:
    pytest.skip("pyopenvdb not installed", allow_module_level=True)

mn.unregister()
mn.register()


@pytest.fixture
def density_file():
    file = data_dir / "emd_24805.map.gz"
    vdb_file = data_dir / "emd_24805.vdb"
    vdb_file.unlink(missing_ok=True)
    # Make all densities are removed
    for o in bpy.data.objects:
        if o.mn.molecule_type == "density":
            bpy.data.objects.remove(o, do_unlink=True)
    return file


def test_density_load(density_file):

    obj = mn.io.density.load(density_file).object
    evaluated = mn.blender.obj.evaluate_using_mesh(obj)
    pos = mn.blender.obj.get_attribute(evaluated, "position")

    assert len(pos) > 1000

    avg = np.mean(pos, axis=0)
    assert np.linalg.norm(avg) > 0.5

    assert obj.mn.molecule_type == "density"
    assert obj.users_collection[0] == mn.blender.coll.mn()


def test_density_centered(density_file):

    # First load using standar parameters to test recreation of vdb
    o = mn.io.density.load(density_file).object
    # Then refresh the scene
    bpy.data.objects.remove(o, do_unlink=True)

    obj = mn.io.density.load(density_file, center=True).object
    evaluated = mn.blender.obj.evaluate_using_mesh(obj)

    pos = mn.blender.obj.get_attribute(evaluated, "position")

    assert len(pos) > 1000

    avg = np.mean(pos, axis=0)
    assert np.linalg.norm(avg) < 0.1


def test_density_invert(density_file):

    # First load using standar parameters to test recreation of vdb
    o = mn.io.density.load(density_file).object
    # Then refresh the scene
    bpy.data.objects.remove(o, do_unlink=True)

    obj = mn.io.density.load(density_file, invert=True).object
    style_node = mn.blender.nodes.get_style_node(obj)
    style_node.inputs["Threshold"].default_value = 0.01
    evaluated = mn.blender.obj.evaluate_using_mesh(obj)

    pos = mn.blender.obj.get_attribute(evaluated, "position")
    # At this threshold after inverting we should have a cube the size of the volume
    assert pos[:, 0].max() > 2.0
    assert pos[:, 1].max() > 2.0
    assert pos[:, 2].max() > 2.0


def test_density_multiple_load():
    file = data_dir / "emd_24805.map.gz"
    obj = mn.io.density.load(file).object
    obj2 = mn.io.density.load(file).object

    assert obj.mn.molecule_type == "density"
    assert obj2.mn.molecule_type == "density"
    assert obj.users_collection[0] == mn.blender.coll.mn()
    assert obj2.users_collection[0] == mn.blender.coll.mn()


@pytest.mark.parametrize('name', ['', 'NewDensity'])
def test_density_naming_op(density_file, name):
    bpy.context.scene.MN_import_density_name = name
    bpy.context.scene.MN_import_density = str(density_file)
    bpy.ops.mn.import_density()

    if name == '':
        object_name = 'emd_24805'
    else:
        object_name = name
    object = bpy.data.objects[object_name]
    assert object
    assert object.name == object_name


@pytest.mark.parametrize('name', ['', 'NewDensity'])
def test_density_naming_api(density_file, name):
    object = mn.io.density.load(density_file, name).object
    if name == '':
        object_name = 'emd_24805'
    else:
        object_name = name

    assert object
    assert object.name == object_name


@pytest.mark.parametrize("invert,node_setup,center", list(itertools.product([True, False], repeat=3)))
def test_density_operator(snapshot, density_file, invert, node_setup, center):
    scene = bpy.context.scene
    scene.MN_import_density = str(density_file)
    scene.MN_import_density_invert = invert
    scene.MN_import_node_setup = node_setup
    scene.MN_import_density_center = center
    scene.MN_import_density_name = ""
    bpy.ops.mn.import_density()
    object = bpy.data.objects['emd_24805']
    snapshot.assert_match(
        sample_attribute_to_string(
            mn.blender.obj.evaluate_using_mesh(object), 'position'),
        "invert_nodesetup_center_positions.txt"
    )
