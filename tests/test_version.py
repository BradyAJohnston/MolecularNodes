import pytest
from addon_helper import get_version
import MolecularNodes as mn
mn.pkg.install()

@pytest.fixture
def bpy_module(cache):
    return cache.get("bpy_module", None)


def test_versionID_pass(bpy_module):
    expect_version = (2, 4, 3)
    return_version = get_version(bpy_module)
    assert expect_version == return_version

def test_versionID_fail(bpy_module):
    expect_version = (0, 1, 1)
    return_version = get_version(bpy_module)
    assert not expect_version == return_version

def test_load_rcsb(bpy_module):
    mol, file = mn.load.open_structure_rcsb('4ozs')
    assert 1 == 1