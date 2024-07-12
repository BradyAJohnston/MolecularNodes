import pytest
import molecularnodes as mn

from .constants import data_dir


@pytest.fixture
def filepath():
    return data_dir / "1f2n.bcif"


def test_bcif_init(filepath):
    bcif = mn.entities.BCIF(filepath)
    assert bcif.file_path == filepath


def test_bcif_read(filepath):
    bcif = mn.entities.BCIF(filepath)
    assert bcif.file is not None


def test_bcif_get_structure(filepath):
    bcif = mn.entities.BCIF(filepath)
    structure = bcif.get_structure()
    assert structure is not None
    # assert structure.shape[1] == len(bcif)


def test_bcif_assemblies(filepath):
    bcif = mn.entities.BCIF(filepath)
    assemblies = bcif.assemblies()
    assert assemblies is not None


def test_bcif_entity_ids(filepath):
    bcif = mn.entities.BCIF(filepath)
    entity_ids = bcif.entity_ids
    assert entity_ids is not None
    assert entity_ids == ["CAPSID PROTEIN", "CALCIUM ION", "water"]
