import numpy as np
import pytest
import molecularnodes as mn

pytest.importorskip("rdkit")


def test_from_smiles():
    mol = mn.Molecule.from_smiles("CCO")
    # ethanol with hydrogens added
    assert len(mol) == 9
    assert sorted(mol.universe.atoms.elements) == [
        "C",
        "C",
        "H",
        "H",
        "H",
        "H",
        "H",
        "H",
        "O",
    ]
    assert len(mol.universe.atoms.bonds) == 8
    assert np.isfinite(mol.universe.atoms.positions).all()
    assert set(mol.named_attribute("atomic_number")) == {1, 6, 8}
    assert mol.name == "CCO"


def test_from_smiles_name_and_style():
    mol = mn.Molecule.from_smiles("c1ccccc1", name="benzene")
    assert mol.name == "benzene"
    mol.add_style("ball_and_stick")


def test_from_smiles_conformers_become_frames():
    mol = mn.Molecule.from_smiles("CCO", numConfs=3)
    assert mol.universe.trajectory.n_frames == 3


def test_from_smiles_invalid_string():
    with pytest.raises(SyntaxError, match="SMILES"):
        mn.Molecule.from_smiles("not_a_smiles")
