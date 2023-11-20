import numpy as np
import pytest
import molecularnodes as mn

def test_read_topology(tmp_path):
    # Create a temporary file with a sample topology
    filepath = tmp_path / "topology.txt"
    with open(filepath, 'w') as file:
        file.write("metadata\n")
        file.write("1 A 2 3\n")
        file.write("4 C 5 6\n")
        file.write("7 G 8 9\n")
    
    # Call the read_topology function
    topology = mn.dna.read_topology(filepath)
    
    # Define the expected topology
    expected_topology = np.array([
        [1, 30, 2, 3],
        [4, 31, 5, 6],
        [7, 32, 8, 9]
    ])
    
    # Assert that the returned topology matches the expected topology
    assert np.array_equal(topology, expected_topology)

def test_topology_to_idx():
    top = np.array([
        [1, 31, -1,  1],
        [1,  3,  0,  1],
        [1,  2,  1, -1]
    ])
    
    bonds = mn.dna.toplogy_to_bond_idx_pairs(top)
    expected = np.array([[0, 1], [2, 1]])
    
    assert np.array_equal(bonds, expected)

def test_base_lookup():
    bases = np.array(['A', 'C', 'C', 'G', 'T', '-10', 'G', 'C', '-3'])
    expected = np.array([30, 31, 31, 32, 33, -1, 32, 31, -1])
    
    ints = mn.dna.base_to_int(bases)
    
    assert np.array_equal(ints, expected)
    