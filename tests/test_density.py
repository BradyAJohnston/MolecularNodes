import molecularnodes as mn
import bpy
import pytest
import numpy as np
from .constants import test_data_directory
from .utils import (
    apply_mods
)
try:
    import pyopenvdb
except ImportError:
    pytest.skip("pyopenvdb not installed", allow_module_level=True)

mn.unregister()
mn.register()

@pytest.fixture
def density_file():
    file = test_data_directory / "emd_24805.map.gz"
    vdb_file = test_data_directory / "emd_24805.vdb"
    vdb_file.unlink(missing_ok=True)
    # Make all densities are removed
    for o in bpy.data.objects:
        if o.mn.molecule_type == "density":
            bpy.data.objects.remove(o,do_unlink=True)
    return file

def test_density_load(density_file):
    
    obj = mn.io.density.load(density_file,style="density_surface")
    evaluated = mn.blender.obj.evaluate_using_mesh(obj)
    pos = mn.blender.obj.get_attribute(evaluated,"position")

    assert len(pos) > 1000

    avg = np.mean(pos,axis=0)
    assert np.linalg.norm(avg) > 0.5 

    assert obj.mn.molecule_type == "density"
    assert obj.users_collection[0] == mn.blender.coll.mn()

def test_density_centered(density_file):

    # First load using standar parameters to test recreation of vdb
    o = mn.io.density.load(density_file,style="density_surface")
    # Then refresh the scene
    bpy.data.objects.remove(o,do_unlink=True)

    obj = mn.io.density.load(density_file,style="density_surface",center=True)
    evaluated = mn.blender.obj.evaluate_using_mesh(obj)

    pos = mn.blender.obj.get_attribute(evaluated,"position")

    assert len(pos) > 1000

    avg = np.mean(pos,axis=0)
    assert np.linalg.norm(avg) < 0.1

def test_density_invert(density_file):

    # First load using standar parameters to test recreation of vdb
    o = mn.io.density.load(density_file,style="density_surface")
    # Then refresh the scene
    bpy.data.objects.remove(o,do_unlink=True)

    obj = mn.io.density.load(density_file,style="density_surface",invert=True)
    style_node = mn.blender.nodes.get_style_node(obj)
    style_node.inputs["Threshold"].default_value = 0.01
    evaluated = mn.blender.obj.evaluate_using_mesh(obj)

    pos = mn.blender.obj.get_attribute(evaluated,"position")    
    # At this threshold after inverting we should have a cube the size of the volume
    assert pos[:,0].max() > 2.0
    assert pos[:,1].max() > 2.0
    assert pos[:,2].max() > 2.0

def test_density_multiple_load():
    file = test_data_directory / "emd_24805.map.gz"
    obj = mn.io.density.load(file,style="density_surface")
    obj2 = mn.io.density.load(file,style="density_surface")

    assert obj.mn.molecule_type == "density"
    assert obj2.mn.molecule_type == "density"
    assert obj.users_collection[0] == mn.blender.coll.mn()
    assert obj2.users_collection[0] == mn.blender.coll.mn()