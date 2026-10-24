import os
import tempfile

# Importing bpy prepends the user's Blender extension wheels (e.g.
# ~/.config/blender/*/extensions/.local/site-packages) to sys.path, shadowing
# this environment's packages with whatever the enabled extensions ship -
# including possibly incomplete copies while an extension is being rebuilt.
# Point bpy at an empty extensions dir so tests only ever see this
# environment's packages. Must be set before bpy is first imported.
os.environ.setdefault(
    "BLENDER_USER_EXTENSIONS", tempfile.mkdtemp(prefix="mn-test-extensions-")
)

import ipaddress  # noqa: E402
import shutil  # noqa: E402
import socket  # noqa: E402
import sys  # noqa: E402
from os.path import dirname, join, realpath  # noqa: E402
from pathlib import Path  # noqa: E402
import bpy  # noqa: E402

# Build the gitignored node-asset library (molecularnodes/assets/nodes.blend)
# when it is missing or the sources under molecularnodes/nodes/ have changed.
# Paths and flags come from [tool.nodebpy.assets] in pyproject.toml; the
# staleness check is a cheap in-process hash, and the build itself needs a
# fresh session so it runs as a subprocess. xdist workers skip it - the
# controller process has already built by the time they spawn.
if os.environ.get("PYTEST_XDIST_WORKER") is None:
    import subprocess  # noqa: E402
    from argparse import Namespace  # noqa: E402
    from nodebpy.assets import _pipeline as _assets  # noqa: E402

    # the CLI parser pre-seeds the positional dests; mirror that here
    _args = Namespace(command="ensure", source=None, blend=None)
    _assets.apply_config(_args, start=Path(__file__).parent)
    _assets.require_positionals(_args)
    if _assets.is_stale(
        _args.blend, _args.source, _args.resources, _assets.stamp_options(_args)
    ):
        _env = os.environ.copy()
        if bpy.app.binary_path:  # inside a full Blender: spawn another
            import nodebpy.assets.__main__ as _assets_main

            # an isolated extensions dir so the spawn can't pick up stale wheels
            _env["BLENDER_USER_EXTENSIONS"] = tempfile.mkdtemp(prefix="mn-assets-ext-")
            _cmd = [
                bpy.app.binary_path,
                "-b",
                "--factory-startup",
                "-P",
                _assets_main.__file__,
                "--",
                "ensure",
            ]
        else:
            _cmd = [sys.executable, "-m", "nodebpy.assets", "ensure"]
        subprocess.run(_cmd, check=True, env=_env, cwd=Path(__file__).parent.parent)
import numpy as np  # noqa: E402
import pytest  # noqa: E402
import molecularnodes as mn  # noqa: E402
from .utils import NumpySnapshotExtension  # noqa: E402

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


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "network: downloads from a remote server; CI runs these in a single job, "
        "deselect with -m 'not network'",
    )


def _is_local(address) -> bool:
    "Whether a socket address is a Unix socket or a loopback address (e.g. IMD tests)."
    if not isinstance(address, tuple):
        return True
    host = address[0]
    if host == "localhost":
        return True
    try:
        return ipaddress.ip_address(host).is_loopback
    except ValueError:
        return False


# Whether the running test may open remote connections; see `no_network`.
_network_allowed = False
_connect = socket.socket.connect


def _guarded_connect(sock, address):
    if not _network_allowed and not _is_local(address):
        raise RuntimeError(
            f"Tried to connect to {address[0]}. Load structures with the `fetch` "
            "fixture (add missing files to tests/data), or mark the test with "
            "@pytest.mark.network."
        )
    return _connect(sock, address)


# patched for the whole session, so module and session scoped fixtures, which are set
# up before any function scoped fixture, can't download either
socket.socket.connect = _guarded_connect


@pytest.fixture(autouse=True)
def no_network(request):
    """
    Refuse remote connections unless the test is marked `network`.

    Structures come from tests/data through the `fetch` fixture, so a test that tries to
    download has a missing data file. Many CI jobs downloading at once get rate limited
    by RCSB and EMDB and fail.
    """
    global _network_allowed
    _network_allowed = request.node.get_closest_marker("network") is not None
    yield
    _network_allowed = False


@pytest.fixture
def fetch():
    """
    Fetch a structure from tests/data, never from the network.

    Takes the same arguments as `mn.Molecule.fetch`. The file must already be in
    tests/data under the code as spelled in the test; download it there once and commit it.
    """

    def _fetch(code: str, format: str = "bcif", **kwargs) -> mn.Molecule:
        format = format.strip(".")
        if not Path(DATA_DIR, f"{code}.{format}").exists():
            raise FileNotFoundError(
                f"{code}.{format} is not in tests/data. Download it there once and "
                "commit it, so tests don't download it on every CI job."
            )
        return mn.Molecule.fetch(code, format=format, cache=DATA_DIR, **kwargs)

    return _fetch


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
