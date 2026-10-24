import os
import tempfile

# see tests/conftest.py, this has to be set before bpy is first imported
os.environ.setdefault(
    "BLENDER_USER_EXTENSIONS", tempfile.mkdtemp(prefix="mn-benchmark-extensions-")
)

import subprocess  # noqa: E402
import sys  # noqa: E402
from argparse import Namespace  # noqa: E402
from pathlib import Path  # noqa: E402
import bpy  # noqa: E402
import pytest  # noqa: E402

# build the gitignored node-asset library if it is missing or stale, as in
# tests/conftest.py
from nodebpy.assets import _pipeline as _assets  # noqa: E402

_args = Namespace(command="ensure", source=None, blend=None)
_assets.apply_config(_args, start=Path(__file__).parent)
_assets.require_positionals(_args)
if _assets.is_stale(
    _args.blend, _args.source, _args.resources, _assets.stamp_options(_args)
):
    subprocess.run(
        [sys.executable, "-m", "nodebpy.assets", "ensure"],
        check=True,
        cwd=Path(__file__).parent.parent,
    )

import molecularnodes as mn  # noqa: E402

mn.ui.addon._test_register()


def reset():
    bpy.ops.wm.read_homefile(app_template="")
    mn.session.get_session().clear()


@pytest.fixture(autouse=True)
def fresh_file():
    reset()
    yield
    reset()


@pytest.fixture
def run(benchmark):
    """
    Benchmark a function, starting each round from a fresh file so data-blocks
    (objects, appended node groups, materials) created in previous rounds don't
    affect the timings.
    """

    def _run(func, *args, rounds: int = 5, **kwargs):
        def setup():
            reset()
            return args, kwargs

        return benchmark.pedantic(func, setup=setup, rounds=rounds)

    return _run
