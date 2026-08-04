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
BLEND_DIR = Path(dirname(realpath(__file__))) / "blend_files"
IS_GITHUB_ACTIONS = os.getenv("GITHUB_ACTIONS") == "true"
IS_SELF_HOSTED = os.getenv("environment") == "self-hosted"


def save_blend_file(request):
    """
    Save the current scene for inspection, named after the test that produced it.

    The module is included in the name so that same-named tests in different files
    don't race on one file under pytest-xdist.
    """
    name = f"{request.module.__name__}.{request.node.name}"
    for char in '/\\:*?"<>|':
        name = name.replace(char, "_")
    BLEND_DIR.mkdir(exist_ok=True)
    # copy=True so the test session's own file state is untouched, and no relative
    # remapping so paths to external files (.vdb, .pdb, .xtc) still resolve
    bpy.ops.wm.save_as_mainfile(
        filepath=str(BLEND_DIR / f"{name}.blend"),
        copy=True,
        relative_remap=False,
    )


@pytest.fixture(autouse=True)
def run_around_tests(request):
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
    # save the scene before it is reset, so a failing test can be opened and inspected
    save_blend_file(request)
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
