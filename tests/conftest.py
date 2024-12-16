import bpy
from os.path import join, dirname, realpath
import sys
import pytest
from .utils import NumpySnapshotExtension
import molecularnodes as mn

# mn.unregister()
# mn.register()


DATA_DIR = join(dirname(realpath(__file__)), "data")


@pytest.fixture(autouse=True)
def run_around_tests():
    # Code that will run before your test, for example:
    bpy.ops.wm.read_homefile(app_template="")
    bpy.context.scene.MNSession.clear()
    # A test function will be run at this point
    assert True
    # Code that will run after your test, for example:
    # files_after = # ... do something to check the existing files
    # assert files_before == files_after


def pytest_sessionstart(session):
    """
    Insert ``MolecularNodes`` into ``PYTHONPATH`` to make it importable.
    """
    project_dir = dirname(dirname(realpath(__file__)))
    sys.path.insert(0, join(project_dir))


@pytest.fixture
def snapshot_custom(snapshot):
    return snapshot.use_extension(NumpySnapshotExtension)
