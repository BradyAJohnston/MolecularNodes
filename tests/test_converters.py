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


# --- Stage 3: the unified loader (structure file -> Universe-backed entity) ---------

# attributes that the unified entity should compute for a file-loaded structure
_EXPECTED_ATTRS = [
    "atomic_number",
    "vdw_radii",
    "mass",
    "res_id",
    "chain_id",
    "atom_name",
    "Color",
    "is_backbone",
    "is_side_chain",
    "is_solvent",
]


@pytest.mark.parametrize("suffix", ["pdb", "cif", "bcif"])
def test_from_file_loads_universe_backed_entity(suffix):
    """A structure file loads into the Universe-backed entity with the full attr set."""
    traj = mn.Trajectory.from_file(data_dir / f"1BNA.{suffix}")

    # it is genuinely Universe-backed
    assert traj.universe.atoms.n_atoms > 0
    assert traj.universe.trajectory.n_frames == 1

    attributes = traj.list_attributes()
    for attr in _EXPECTED_ATTRS:
        assert attr in attributes, f"missing attribute '{attr}'"

    # metadata parsed from the file is stored on the object
    assert traj.props.chain_ids == ["A", "B"]
    assert traj.props.filepath.endswith(f"1BNA.{suffix}")


def test_from_file_selections_work():
    """MDAnalysis selection strings work on a file-loaded structure."""
    traj = mn.Trajectory.from_file(data_dir / "1BNA.pdb")
    item = traj.selections.from_string("resid 1-4")
    assert item.name in traj.list_attributes()
    # the boolean selection actually selects a subset of atoms
    mask = traj.named_attribute(item.name)
    assert 0 < mask.sum() < traj.universe.atoms.n_atoms


def test_from_file_sdf_has_bonds():
    """A small-molecule SDF loads with connectivity (edges) intact."""
    traj = mn.Trajectory.from_file(data_dir / "caffeine.sdf")
    assert len(traj.data.edges) > 0


@pytest.mark.parametrize("code", multimodel_codes)
def test_from_file_multimodel_frames(code):
    """A multi-model structure file becomes a multi-frame Universe-backed entity."""
    n_models = _read_array(code).stack_depth()
    traj = mn.Trajectory.from_file(data_dir / f"{code}.bcif")
    assert traj.universe.trajectory.n_frames == n_models


def test_fetch_stores_source_and_assemblies():
    """fetch records the code/database and exposes biological assemblies."""
    traj = mn.Trajectory.fetch("4ozs", format=".bcif", cache=data_dir)
    assert traj.props.code == "4ozs"
    assert traj.props.database == "rcsb"
    assemblies = traj.assemblies()
    assert assemblies is not None and len(assemblies) > 0
