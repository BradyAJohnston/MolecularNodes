import os
import shutil
import sys
from os.path import dirname, join, realpath
from pathlib import Path
import bpy
import numpy as np
import pytest
import molecularnodes as mn
from .utils import NumpySnapshotExtension

# Pin numpy print format so snapshots are consistent across numpy 1.x and 2.x
# (numpy 2.2+ adds shape= to array_repr for truncated arrays)
# TODO: Remove legacy="1.25" when bpy 5.1.1 is released
if int(np.__version__.split(".")[0]) >= 2:
    np.set_printoptions(legacy="2.1")

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


@pytest.fixture
def isolated_density_file(tmp_path):
    """
    Copy a density file into a per-test temporary directory.

    The generated `.vdb` is written next to the file it was created from, so tests
    sharing a file in `tests/data` overwrite and delete each other's `.vdb` when run
    in parallel. Copying gives each test its own directory to write into.
    """

    def _copy(file: str | Path) -> Path:
        file = Path(file)
        destination = tmp_path / file.name
        shutil.copy(file, destination)
        return destination

    return _copy
