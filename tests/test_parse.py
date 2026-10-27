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


def test_entity_ids_with_commas():
    # 6ZPA has an entity described as "(2R,3S,5R)-5-(6-amino-...", which a
    # comma-joined string split into several entity IDs
    from molecularnodes.entities.molecule.reader import read_structure

    path = data_dir / "6ZPA.bcif"
    mol = mn.Molecule.load(path)
    assert mol.props.entity_ids == read_structure(path).entity_ids()
    assert any("," in entity_id for entity_id in mol.props.entity_ids)


def test_comma_joined_ids_still_read(filepath):
    # files saved before these were stored as JSON hold comma-joined strings
    mol = mn.Molecule.load(filepath)
    mol.props.internal_chain_ids = "A,B,C"
    assert mol.props.chain_ids == ["A", "B", "C"]
