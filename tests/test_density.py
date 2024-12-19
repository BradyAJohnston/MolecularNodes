import molecularnodes as mn
import bpy
import pytest
import numpy as np
import itertools
from .constants import data_dir
from .utils import NumpySnapshotExtension
from databpy import ObjectTracker

try:
    import pyopenvdb
except ImportError:
    pytest.skip("pyopenvdb not installed", allow_module_level=True)


@pytest.fixture
def density_file():
    file = data_dir / "emd_24805.map.gz"
    vdb_file = data_dir / "emd_24805.vdb"
    vdb_file.unlink(missing_ok=True)
    # Make all densities are removed
    for o in bpy.data.objects:
        if o.mn.entity_type == "density":
            bpy.data.objects.remove(o, do_unlink=True)
    return file


def test_density_load(density_file):
    density = mn.entities.density.load(density_file)
    pos = density.named_attribute("position")

    assert len(pos) > 1000

    avg = np.mean(pos, axis=0)
    assert np.linalg.norm(avg) > 0.5

    print(f"{list(bpy.data.objects)=}")

    assert density.object.mn.entity_type == "density"
    assert density.object.users_collection[0] == mn.blender.coll.mn()


def test_density_centered(density_file):
    # First load using standard parameters to test recreation of vdb
    density = mn.entities.density.load(density_file, center=True, overwrite=True)
    pos = density.named_attribute("position")
    avg = np.mean(pos, axis=0)
    assert len(pos) > 1000
    assert np.linalg.norm(avg) < 0.1


def test_density_invert(density_file):
    # First load using standard parameters to test recreation of vdb
    density = mn.entities.density.load(density_file)

    # Then refresh the scene
    bpy.data.objects.remove(density.object, do_unlink=True)

    density = mn.entities.density.load(density_file, invert=True)
    style_node = mn.blender.nodes.get_style_node(density.object)
    style_node.inputs["Threshold"].default_value = 0.01

    pos = density.named_attribute("position")
    # At this threshold after inverting we should have a cube the size of the volume
    assert pos[:, 0].max() > 2.0
    assert pos[:, 1].max() > 2.0
    assert pos[:, 2].max() > 2.0


def test_density_multiple_load():
    file = data_dir / "emd_24805.map.gz"
    density1 = mn.entities.density.load(file)
    density2 = mn.entities.density.load(file)

    assert density1.object.mn.entity_type == "density"
    assert density2.object.mn.entity_type == "density"
    assert density1.object.users_collection[0] == mn.blender.coll.mn()
    assert density2.object.users_collection[0] == mn.blender.coll.mn()


@pytest.mark.parametrize("name", ["", "NewDensity"])
def test_density_naming_op(density_file, name):
    bpy.context.scene.MN_import_density_name = name
    bpy.context.scene.MN_import_density = str(density_file)
    bpy.ops.mn.import_density()

    if name == "":
        object_name = "emd_24805"
    else:
        object_name = name
    object = bpy.data.objects[object_name]
    assert object
    assert object.name == object_name


@pytest.mark.parametrize("name", ["", "NewDensity"])
def test_density_naming_api(density_file, name):
    density = mn.entities.density.load(density_file, name)

    if name == "":
        object_name = "emd_24805"
    else:
        object_name = name

    assert density.object
    assert density.object.name == object_name


@pytest.mark.parametrize(
    "invert,node_setup,center", list(itertools.product([True, False], repeat=3))
)
def test_density_operator(
    snapshot_custom: NumpySnapshotExtension, density_file, invert, node_setup, center
):
    scene = bpy.context.scene
    scene.MN_import_density = str(density_file)
    scene.MN_import_density_invert = invert
    scene.MN_import_node_setup = node_setup
    scene.MN_import_density_center = center
    scene.MN_import_density_name = ""
    with ObjectTracker() as o:
        bpy.ops.mn.import_density()
        density = scene.MNSession.match(o.latest())

    assert snapshot_custom == density.position
