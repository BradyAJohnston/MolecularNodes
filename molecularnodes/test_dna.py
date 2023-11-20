import numpy as np
import pytest
from molecularnodes.dna import read_topology, toplogy_to_bond_idx_pairs

def test_read_topology(tmp_path):
    # Create a temporary file with a sample topology
    filepath = tmp_path / "topology.txt"
    with open(filepath, 'w') as file:
        file.write("metadata\n")
        file.write("1 A 2 3\n")
        file.write("4 C 5 6\n")
        file.write("7 G 8 9\n")
    
    # Call the read_topology function
    topology = read_topology(filepath)
    
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
    
    bonds = toplogy_to_bond_idx_pairs(top)
    expected = np.array([[0, 1], [2, 1]])
    
    assert np.array_equal(bonds, expected)