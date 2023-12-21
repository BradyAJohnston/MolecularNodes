import molecularnodes as mn
import pytest
import numpy as np
from .constants import test_data_directory
try:
    import pyopenvdb
except ImportError:
    pytest.skip("pyopenvdb not installed", allow_module_level=True)

mn.unregister()
mn.register()

def test_density_load():
    file = test_data_directory / "emd_24805.map.gz"
    obj = mn.io.density.load(file,style="density_surface")

    assert obj.mn.molecule_type == "density"
    assert obj.users_collection[0] == mn.blender.coll.mn()

def test_density_multiple_load():
    file = test_data_directory / "emd_24805.map.gz"
    obj = mn.io.density.load(file,style="density_surface")
    obj2 = mn.io.density.load(file,style="density_surface")

    assert obj.mn.molecule_type == "density"
    assert obj2.mn.molecule_type == "density"
    assert obj.users_collection[0] == mn.blender.coll.mn()
    assert obj2.users_collection[0] == mn.blender.coll.mn()