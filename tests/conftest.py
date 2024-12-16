import bpy
from os.path import join, dirname, realpath
import sys
import pytest
from .utils import NumpySnapshotExtension
import molecularnodes as mn


DATA_DIR = join(dirname(realpath(__file__)), "data")


@pytest.fixture(autouse=True)
def run_around_tests():
    # Code that will run before each tests

    bpy.ops.wm.read_homefile(app_template="")
    mn.session.get_session().clear()
    for tree in bpy.data.node_groups:
        bpy.data.node_groups.remove(tree)

    print(f"{mn.session.get_session().entities=}")
    print(f"{list(bpy.data.objects)=}")

    yield

    bpy.ops.wm.read_homefile(app_template="")
    mn.session.get_session().clear()
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
