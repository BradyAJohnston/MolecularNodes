import itertools
import MDAnalysis as mda
import numpy as np
import pytest
import molecularnodes as mn
from .constants import codes, data_dir


@pytest.mark.parametrize(
    "code, format", itertools.product(codes, ["bcif", "cif", "pdb"])
)
def test_biotite_converter(code, format):
    # fetch and load molecule through biotite
    mol = mn.Molecule.fetch(code, format=format, cache=data_dir)
    # create MDAnalysis universe from biotite structure
    u = mda.Universe(mn.converters.BiotiteWrapper(mol.array))
    # load trajectory into Blender
    traj = mn.Trajectory(u)
    # check coords
    assert np.allclose(mol.position, traj.position)
    # check named attributes
    attrs = [
        "chain_id",
        "res_id",
        "res_name",
        "atom_name",
        "atom_id",
        "b_factor",
        "occupancy",
        "charge",
    ]
    for attr in attrs:
        assert np.allclose(mol.named_attribute(attr), traj.named_attribute(attr))
    # check string attributes
    assert np.array_equal(mol.array.element, traj.elements)
    assert np.array_equal(mol.array.ins_code, traj.atoms.icodes)
    # check computed attributes
    # the ones that currently differ are commented out with reasons
    computed_attrs = [
        ("mass", 1e-3),
        # ("atomic_number", 0),  # unsure which is correct, computation different
        # ("vdw_radii", 0), # np.char.title in Molecule causes difference
        ("charge", 0),
        ("is_alpha_carbon", 0),
        ("is_solvent", 0),
        # ("is_backbone", 0),  # computation different
        ("is_nucleic", 0),
        ("is_peptide", 0),
        # ("ures_id", 0),  # not currently in Trajectory
        # ("lipophobicity", 0),  # not currently in Trajectory
        # ("Color", 0),  # not currently in Trajectory
        # ("is_hetero", 0),  # not currently in Trajectory
        # ("is_side_chain", 0),  # not currently in Trajectory
        # ("is_carb", 0),  # not currently in Trajectory
        # ("sec_struct", 0),  # not currently in Trajectory
        # ("entity_id", 0),  # not currently in Trajectory
        # ("asym_id", 0),  # not currently in Trajectory
        # ("pdb_model_num", 0),  # not currently in Trajectory
    ]
    for attr, rtol in computed_attrs:
        assert np.allclose(
            mol.named_attribute(attr), traj.named_attribute(attr), rtol=rtol
        )
