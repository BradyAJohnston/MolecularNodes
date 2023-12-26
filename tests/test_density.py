import molecularnodes as mn
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

def test_density_load():
    file = test_data_directory / "emd_24805.map.gz"
    vdb_file = test_data_directory / "emd_24805.vdb"
    vdb_file.unlink(missing_ok=True)
    obj = mn.io.density.load(file,style="density_surface")
    evaluated = mn.blender.obj.evaluate_using_debug_cube(obj)

    pos = mn.blender.obj.get_attribute(evaluated,"position")

    assert len(pos) > 1000
    # Get the average of positions
    avg = np.mean(pos,axis=0)
    
    assert np.linalg.norm(avg) > 0.1    
    assert obj.mn.molecule_type == "density"
    assert obj.users_collection[0] == mn.blender.coll.mn()

def test_density_centered():
    file = test_data_directory / "emd_24805.map.gz"
    vdb_file = test_data_directory / "emd_24805.vdb"
    vdb_file.unlink(missing_ok=True)
    obj = mn.io.density.load(file,style="density_surface")
    # Need to load an empty file to reset vdb file
    mn.bpy.ops.wm.read_factory_settings()
    obj2 = mn.io.density.load(file,style="density_surface",center=True)
    evaluated = mn.blender.obj.evaluate_using_debug_cube(obj2)

    pos = mn.blender.obj.get_attribute(evaluated,"position")

    assert len(pos) > 1000
    # Get the average of positions
    avg = np.mean(pos,axis=0)
    
    assert np.linalg.norm(avg) < 0.1
    assert not np.allclose(avg,[0.0,0.0,0.0])

def test_density_multiple_load():
    file = test_data_directory / "emd_24805.map.gz"
    obj = mn.io.density.load(file,style="density_surface")
    obj2 = mn.io.density.load(file,style="density_surface")

    assert obj.mn.molecule_type == "density"
    assert obj2.mn.molecule_type == "density"
    assert obj.users_collection[0] == mn.blender.coll.mn()
    assert obj2.users_collection[0] == mn.blender.coll.mn()