import pytest
import molecularnodes as mn

from .constants import data_dir


@pytest.fixture
def filepath():
    return data_dir / '1f2n.mmtf'


def test_mmtf_init(filepath):
    mmtf = mn.io.parse.MMTF(filepath)
    assert mmtf.file_path == filepath


def test_mmtf_read(filepath):
    mmtf = mn.io.parse.MMTF(filepath)
    assert mmtf.file is not None


def test_mmtf_get_structure(filepath):
    mmtf = mn.io.parse.MMTF(filepath)
    structure = mmtf._get_structure()
    assert structure is not None
    assert structure.shape[0] == mmtf.n_models
    assert structure.shape[1] == mmtf.n_atoms


def test_mmtf_assemblies(filepath):
    mmtf = mn.io.parse.MMTF(filepath)
    assemblies = mmtf._assemblies()
    assert assemblies is not None


def test_mmtf_entity_ids(filepath):
    mmtf = mn.io.parse.MMTF(filepath)
    entity_ids = mmtf.entity_ids
    assert entity_ids is not None
    assert entity_ids == ['CAPSID PROTEIN', 'CALCIUM ION', 'water']
