import itertools
import bpy
import numpy as np
import pytest
from databpy import ObjectTracker
import molecularnodes as mn
from molecularnodes.nodes import nodes
from .constants import data_dir
from .utils import NumpySnapshotExtension


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
    style_node = nodes.get_style_node(density.object)
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


def test_density_naming_op(density_file):
    bpy.context.scene.mn.import_density = str(density_file)
    bpy.ops.mn.import_density()
    _object = bpy.data.objects[density_file.name]


@pytest.mark.parametrize(
    "invert,node_setup,center", list(itertools.product([True, False], repeat=3))
)
def test_density_operator(
    snapshot_custom: NumpySnapshotExtension, density_file, invert, node_setup, center
):
    scene = bpy.context.scene
    scene.mn.import_density = str(density_file)
    scene.mn.import_density_invert = invert
    scene.mn.import_node_setup = node_setup
    scene.mn.import_density_center = center
    with ObjectTracker() as o:
        bpy.ops.mn.import_density()
        density: mn.entities.Density = scene.MNSession.match(o.latest())

    obj = mn.blender.mesh.evaluate_using_mesh(density.object)
    assert len(obj.data.vertices) == [198862, 287422][int(invert)]


@pytest.fixture
def density_file_dx():
    file = data_dir / "water.dx.gz"
    vdb_file = data_dir / "water.vdb"
    vdb_file.unlink(missing_ok=True)
    # Make sure all densities are removed
    for o in bpy.data.objects:
        if o.mn.entity_type == "density":
            bpy.data.objects.remove(o, do_unlink=True)
    return file


def test_density_load_dx(density_file_dx):
    density = mn.entities.density.load(density_file_dx)
    print(f"{list(bpy.data.objects)=}")
    assert density.object.mn.entity_type == "density"
    assert density.object.users_collection[0] == mn.blender.coll.mn()
