import numpy as np
import pytest
from biotite import InvalidFileError
import molecularnodes as mn
from molecularnodes.entities.molecule import corecif
from .constants import data_dir

CRYSTAL_FILE = data_dir / "mn_test_crystal.cif"
SF_FILE = data_dir / "structure_factors.cif"


def test_is_core_cif():
    assert corecif.is_core_cif(CRYSTAL_FILE)
    # mmCIF and structure factor files use dotted category.column tags
    assert not corecif.is_core_cif(data_dir / "8H1B.cif")
    assert not corecif.is_core_cif(SF_FILE)


def test_parse_core_cif():
    crystal = corecif.parse_core_cif(CRYSTAL_FILE)

    # standard uncertainty suffixes are dropped, including the empty '()'
    # some converters write (cell 4.000(2), Cl1 fract_z 0.25())
    assert np.allclose(crystal.cell, [4.0, 5.0, 6.0, 90.0, 90.0, 90.0])

    # Na on the origin maps onto itself under inversion (2 unique positions);
    # Cl on a general position expands under all four operators
    assert list(crystal.names) == ["Na1"] * 2 + ["Cl1"] * 4
    assert list(crystal.elements) == ["Na"] * 2 + ["Cl"] * 4
    # missing occupancy ('.') defaults to fully occupied
    assert np.allclose(crystal.occupancies, 1.0)

    expected_fract = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.25, 0.25, 0.25],
            [0.25, 0.25, 0.75],
            [0.75, 0.75, 0.25],
            [0.75, 0.75, 0.75],
        ]
    )
    expected = expected_fract * crystal.cell[:3]
    assert np.allclose(
        sorted(map(tuple, crystal.positions)), sorted(map(tuple, expected)), atol=1e-4
    )


def test_load_core_cif():
    mol = mn.Molecule.load(CRYSTAL_FILE)
    assert len(mol) == 6
    assert set(mol.universe.atoms.elements) == {"Na", "Cl"}
    assert set(mol.named_attribute("atomic_number")) == {11, 17}


def test_structure_factor_file_errors():
    with pytest.raises(InvalidFileError, match="structure factor"):
        mn.Molecule.load(SF_FILE)


def test_load_xyz():
    mol = mn.Molecule.load(data_dir / "mn_test.xyz")
    assert len(mol) == 5
    assert list(mol.universe.atoms.elements) == ["C", "H", "H", "H", "H"]
    assert set(mol.named_attribute("atomic_number")) == {6, 1}
