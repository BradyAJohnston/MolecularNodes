import pytest
print("place 9")
from addon_helper import get_version
print("place 10")
import MolecularNodes as mn
print("place 11")



@pytest.fixture
def bpy_module(cache):
    return cache.get("bpy_module", None)

print("place 12")
# ensure we can successfully install all of the required pacakges 
def test_install_packages(bpy_module):
    mn.pkg.install_all_packages()
    assert 1 == 1

print("place 13")
def test_versionID_pass(bpy_module):
    expect_version = (2, 5, 2)
    return_version = get_version(bpy_module)
    assert expect_version == return_version
print("place 14")

def test_versionID_fail(bpy_module):
    expect_version = (0, 1, 1)
    return_version = get_version(bpy_module)
    assert not expect_version == return_version

# able to load a structure without any trouble
def test_load_rcsb(bpy_module):
    mol, file = mn.load.open_structure_rcsb('4ozs')
    assert 1 == 1