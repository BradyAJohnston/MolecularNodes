import itertools
import numpy as np
import pytest
import molecularnodes as mn
from molecularnodes.converters import universe_from_atoms
from molecularnodes.entities.molecule.base import Molecule
from .constants import codes, data_dir

# multi-model (NMR ensemble) structures cached under tests/data
multimodel_codes = ["1NMR", "2M6Q"]


def _read_array(code):
    """Read a cached structure into a biotite AtomArrayStack (no Blender entity)."""
    return Molecule._read(data_dir / f"{code}.bcif").array


@pytest.mark.parametrize(
    "code, format", list(itertools.product(codes, ["bcif", "cif", "pdb"]))
)
def test_biotite_converter(code, format):
    # fetch and load molecule through biotite
    mol = mn.Molecule.fetch(code, format=format, cache=data_dir)
    # create MDAnalysis universe from biotite structure (via the full converter, which
    # carries bonds and file-parsed extras across)
    u = universe_from_atoms(mol.array)
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
    assert np.array_equal(mol.array.element, traj.atoms.elements)
    assert np.array_equal(mol.array.ins_code, traj.atoms.icodes)
    # check computed attributes: these are now recomputed MDAnalysis-side (or carried
    # across the converter) so that both backends produce matching named attributes
    computed_attrs = [
        ("mass", 1e-3),
        ("atomic_number", 0),
        ("vdw_radii", 0),
        ("charge", 0),
        ("is_alpha_carbon", 0),
        ("is_solvent", 0),
        ("is_backbone", 0),
        ("is_nucleic", 0),
        ("is_peptide", 0),
        ("ures_id", 0),
        ("lipophobicity", 0),
        ("Color", 0),
        ("is_hetero", 0),
        ("is_side_chain", 0),
        ("is_carb", 0),
        ("sec_struct", 0),
        ("entity_id", 0),
        # ("asym_id", 0),  # not parsed into the biotite array; revisit if needed
        # ("pdb_model_num", 0),  # not parsed into the biotite array; revisit if needed
    ]
    # some attributes are structure-dependent (e.g. sec_struct only exists for proteins);
    # only compare the ones the biotite backend actually produced for this structure
    mol_attrs = mol.list_attributes(drop_hidden=False)
    for attr, rtol in computed_attrs:
        if attr not in mol_attrs:
            continue
        assert np.allclose(
            mol.named_attribute(attr), traj.named_attribute(attr), rtol=rtol
        ), f"attribute '{attr}' differs between biotite and MDAnalysis backends"


@pytest.mark.parametrize("code", codes)
def test_universe_from_atoms_topology(code):
    """The converter faithfully carries topology, bonds and file-parsed extras."""
    stack = _read_array(code)
    array = stack[0]
    u = universe_from_atoms(stack)

    # single-model structures become a single-frame universe
    assert u.trajectory.n_frames == 1
    assert np.allclose(u.atoms.positions, array.coord)

    # standard topology attributes match the source array
    assert np.array_equal(u.atoms.elements, array.element)
    assert np.array_equal(u.atoms.icodes, array.ins_code)
    assert np.array_equal(u.atoms.resids, array.res_id)
    assert np.array_equal(u.atoms.resnames, array.res_name)
    assert np.array_equal(u.atoms.names, array.atom_name)

    # bonds are carried across so connectivity-based styles work
    assert array.bonds is not None
    assert len(u.atoms.bonds) == array.bonds.as_array().shape[0]

    # file-parsed, non-recomputable annotations ride across as custom attributes
    for name, attr in (
        ("entity_id", "entity_ids"),
        ("sec_struct", "sec_structs"),
        ("hetero", "heteros"),
    ):
        if name in array.get_annotation_categories():
            assert np.array_equal(getattr(u.atoms, attr), array.get_annotation(name))


@pytest.mark.parametrize("code", multimodel_codes)
def test_universe_from_atoms_multimodel(code):
    """A multi-model stack becomes a multi-frame universe with matching coordinates."""
    stack = _read_array(code)
    n_models = stack.stack_depth()
    assert n_models > 1  # guard: these fixtures really are ensembles

    u = universe_from_atoms(stack)
    assert u.trajectory.n_frames == n_models

    # every trajectory frame matches the corresponding model's coordinates
    for i, _ in enumerate(u.trajectory):
        assert np.allclose(u.atoms.positions, stack.coord[i])

    # topology (shared across models) is still intact, including bonds
    assert np.array_equal(u.atoms.elements, stack[0].element)
    if stack[0].bonds is not None:
        assert len(u.atoms.bonds) == stack[0].bonds.as_array().shape[0]
