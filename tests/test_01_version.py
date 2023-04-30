from addon_helper import get_version
import MolecularNodes as mn
import pytest

@pytest.fixture
def bpy_module(cache):
    return cache.get("bpy_module", None)

# ensure we can successfully install all of the required pacakges 
def test_install_packages(bpy_module):
    mn.pkg.install_all_packages()
    assert mn.pkg.is_current('biotite') == True

def test_versionID_pass(bpy_module):
    expect_version = (2, 6, 2)
    return_version = get_version(bpy_module)
    assert expect_version == return_version

def test_versionID_fail(bpy_module):
    expect_version = (2, 5, 0)
    return_version = get_version(bpy_module)
    assert not expect_version == return_version