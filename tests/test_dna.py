import numpy as np
import pytest
import molecularnodes as mn
from molecularnodes.io import dna
from .utils import (
    evaluate,
    sample_attribute_to_string
)
from .constants import (
    data_dir
)


def test_read_topology():
    file_new = data_dir / "oxdna/top_new.top"
    file_old = data_dir / "oxdna/top_old.top"

    arr_old = dna.toplogy_to_bond_idx_pairs(dna.read_topology_old(file_old))
    arr_new = dna.toplogy_to_bond_idx_pairs(dna.read_topology_new(file_new))

    assert np.array_equal(arr_old, arr_new)


def test_topology_to_idx():
    top = np.array([
        [1, 31, -1,  1],
        [1,  3,  0,  1],
        [1,  2,  1, -1]
    ])

    bonds = dna.toplogy_to_bond_idx_pairs(top)
    expected = np.array([[0, 1], [1, 2]])

    assert np.array_equal(bonds, expected)


def test_base_lookup():
    bases = np.array(['A', 'C', 'C', 'G', 'T', '-10', 'G', 'C', '-3'])
    expected = np.array([30, 31, 31, 32, 33, -1, 32, 31, -1])

    ints = dna.base_to_int(bases)

    assert np.array_equal(ints, expected)


def test_read_trajectory():
    traj = dna.read_trajectory(data_dir / "oxdna/holliday.dat")

    assert traj.shape == (20, 98, 15)


def test_read_oxdna(snapshot):
    name = 'holliday'
    mol, coll_frames = dna.load(
        top=data_dir / "oxdna/holliday.top",
        traj=data_dir / "oxdna/holliday.dat",
        name=name
    )

    assert len(coll_frames.objects) == 20
    assert mol.name == name

    for att in mol.data.attributes.keys():
        snapshot.assert_match(
            sample_attribute_to_string(mol, att),
            f"mesh_att_{att}_values.txt"
        )

    # realise all of the geometry and sample some attributes
    mn.blender.nodes.realize_instances(mol)
    for att in mol.data.attributes.keys():
        snapshot.assert_match(
            sample_attribute_to_string(evaluate(mol), att),
            f"realized_mesh_att_{att}_values.txt"
        )
