import pytest
import molecularnodes as mn
from .constants import data_dir


@pytest.fixture
def filepath():
    return data_dir / "1f2n.bcif"


def test_bcif_init(filepath):
    from molecularnodes.entities.molecule.reader import read_structure

    reader = read_structure(filepath)
    assert reader.file_path == filepath


def test_bcif_assemblies(filepath):
    bcif = mn.Molecule.load(filepath)
    assemblies = bcif.assemblies()
    assert assemblies is not None


def test_bcif_entity_ids(filepath):
    mol = mn.Molecule.load(filepath)
    entity_ids = list(mol.props.entity_ids)
    assert entity_ids is not None
    assert entity_ids == ["CAPSID PROTEIN", "CALCIUM ION", "water"]
