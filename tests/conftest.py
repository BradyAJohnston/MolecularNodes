from os.path import join, dirname, realpath
import sys
import pytest
from .utils import NumpySnapshotExtension


DATA_DIR = join(dirname(realpath(__file__)), "data")


def pytest_sessionstart(session):
    """
    Insert ``MolecularNodes`` into ``PYTHONPATH`` to make it importable.
    """
    project_dir = dirname(dirname(realpath(__file__)))
    sys.path.insert(0, join(project_dir))


@pytest.fixture
def snapshot_custom(snapshot):
    return snapshot.use_extension(NumpySnapshotExtension)
