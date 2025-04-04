import os
import sys
from os.path import dirname, join, realpath
import bpy
import pytest
import molecularnodes as mn
from .utils import NumpySnapshotExtension

mn.ui.addon._test_register()


DATA_DIR = join(dirname(realpath(__file__)), "data")
IS_GITHUB_ACTIONS = os.getenv("GITHUB_ACTIONS") == "true"
IS_SELF_HOSTED = os.getenv("environment") == "self-hosted"


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
    print("Post Test setup:")
    print(f"{bpy.app.handlers.frame_change_pre=}")
    print(f"{mn.session.get_session().entities.keys()=}")
    print(f"{list(bpy.data.objects)=}")
    print(f"{list(o.uuid for o in bpy.data.objects)=}")
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
