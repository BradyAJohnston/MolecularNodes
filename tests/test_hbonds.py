import pickle
import bpy
import MDAnalysis as mda
import numpy as np
import pytest
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from nodebpy import geometry as g
import molecularnodes as mn
from molecularnodes.nodes.geometry import FindHBonds, StyleHBonds
from .constants import data_dir
from .utils import mn_ops

pytestmark = [
    pytest.mark.filterwarnings("ignore:.*Empty string to select atoms.*:UserWarning"),
    pytest.mark.filterwarnings(
        "ignore:.*Unknown masses are set to.*:PendingDeprecationWarning"
    ),
]


@pytest.fixture
def universe():
    # gro topology: no bonds or elements, so the name and distance based
    # defaults are used
    return mda.Universe(
        data_dir / "md_ppr/box.gro", data_dir / "md_ppr/first_5_frames.xtc"
    )


def edges(obj: bpy.types.Object | None) -> np.ndarray:
    assert obj is not None
    mesh = obj.data
    assert isinstance(mesh, bpy.types.Mesh)
    pairs = np.zeros(len(mesh.edges) * 2, dtype=int)
    mesh.edges.foreach_get("vertices", pairs)
    return pairs.reshape(-1, 2)


def expected_pairs(universe: mda.Universe, frame: int, **kwargs) -> np.ndarray:
    hba = HydrogenBondAnalysis(universe, update_selections=False, **kwargs)
    hba.run(start=frame, stop=frame + 1)
    return hba.results.hbonds[:, 2:4].astype(int)


def evaluated_vertex_count(obj: bpy.types.Object) -> int:
    depsgraph = bpy.context.evaluated_depsgraph_get()
    mesh = obj.evaluated_get(depsgraph).data
    assert isinstance(mesh, bpy.types.Mesh)
    return len(mesh.vertices)


def test_hbonds_trajectory(universe):
    mol = mn.Molecule(universe)
    obj = mol.calculations.hbonds().object
    assert obj is not None
    assert isinstance(obj.data, bpy.types.Mesh)
    assert mol.frame_attributes["hbonds"] == "CALCULATION"
    assert mol.props.calculations["hbonds"].object == obj
    assert obj.name == f".data_{mol.name}_hbonds"
    assert len(obj.data.vertices) == len(universe.atoms)
    # an independent universe so running the analysis doesn't move the other
    u = mda.Universe(
        data_dir / "md_ppr/box.gro", data_dir / "md_ppr/first_5_frames.xtc"
    )
    for frame in [1, 4, 2]:
        mol.set_frame(frame)
        # the analysis rewinds the trajectory, which must be restored
        assert universe.trajectory.frame == frame
        expected = expected_pairs(
            u,
            frame,
            donors_sel="name N* O*",
            hydrogens_sel="name H*",
            acceptors_sel="name N* O*",
        )
        assert len(expected) > 0
        np.testing.assert_array_equal(edges(obj), expected)


def test_hbonds_structure(fetch):
    # NMR ensemble with hydrogens, bonds and elements
    mol = fetch("2M6Q")
    obj = mol.calculations.hbonds(selection="protein").object
    u = mol.universe.copy()
    for frame in [0, 3]:
        mol.set_frame(frame)
        expected = expected_pairs(
            u,
            frame,
            donors_sel="protein and element N O",
            hydrogens_sel="protein and element H D",
            acceptors_sel="protein and element N O",
        )
        assert len(expected) > 0
        np.testing.assert_array_equal(edges(obj), expected)


def test_hbonds_no_hydrogens(fetch):
    mol = fetch("1l58")
    with pytest.raises(ValueError, match="No hydrogens found"):
        mol.calculations.hbonds()
    assert "hbonds" not in mol.calculations


def test_hbonds_not_updating(fetch):
    mol = fetch("2M6Q")
    obj = mol.calculations.hbonds(updating=False).object
    pairs = edges(obj)
    mol.set_frame(3)
    np.testing.assert_array_equal(edges(obj), pairs)
    # changing a setting recalculates straight away
    mol.calculations.hbonds(d_a_cutoff=3.5)
    assert not np.array_equal(edges(obj), pairs)
    # cutoffs are stored as properties, not as other arguments
    assert mol.calculations._hbonds._kwargs == {}
    assert mol.props.hbonds.d_a_cutoff == pytest.approx(3.5)
    mol.calculations.hbonds(updating=True)
    mol.set_frame(0)
    u = mol.universe.copy()
    np.testing.assert_array_equal(
        edges(obj),
        expected_pairs(
            u,
            0,
            d_a_cutoff=3.5,
            donors_sel="(all) and (element N O)",
            hydrogens_sel="(all) and (element H D)",
            acceptors_sel="(all) and (element N O)",
        ),
    )


def test_hbonds_cutoff_properties(fetch):
    mol = fetch("2M6Q")
    hbonds = mol.calculations.hbonds(d_h_a_angle_cutoff=140)
    props = mol.props.hbonds
    assert np.degrees(props.d_h_a_angle_cutoff) == pytest.approx(140)
    assert hbonds.cutoffs == pytest.approx(
        {"d_a_cutoff": 3.0, "d_h_a_angle_cutoff": 140, "d_h_cutoff": 1.2}
    )
    pairs = edges(hbonds.object)
    # changing a property, as in the UI, recalculates straight away
    props.d_a_cutoff = 3.5
    props.d_h_a_angle_cutoff = np.radians(130)
    new_pairs = edges(hbonds.object)
    assert len(new_pairs) > len(pairs)
    expected = expected_pairs(
        mol.universe.copy(),
        mol.universe.trajectory.frame,
        d_a_cutoff=3.5,
        d_h_a_angle_cutoff=130,
        donors_sel="(all) and (element N O)",
        hydrogens_sel="(all) and (element H D)",
        acceptors_sel="(all) and (element N O)",
    )
    np.testing.assert_array_equal(new_pairs, expected)
    # the selection and other arguments are kept
    mol.calculations.hbonds(selection="protein")
    assert hbonds.cutoffs["d_a_cutoff"] == pytest.approx(3.5)


def test_hbonds_after_session_save(universe, tmp_path):
    # saving the session replaces the universe's trajectory reader, and the
    # analysis reads frames through the new one
    mol = mn.Molecule(universe)
    obj = mol.calculations.hbonds().object
    reader = mol.universe.trajectory
    mn.session.get_session().pickle(tmp_path / "test.blend")
    assert mol.universe.trajectory is not reader
    u = mda.Universe(
        data_dir / "md_ppr/box.gro", data_dir / "md_ppr/first_5_frames.xtc"
    )
    for frame in [2, 3]:
        mol.set_frame(frame)
        expected = expected_pairs(
            u,
            frame,
            donors_sel="(all) and (name N* O*)",
            hydrogens_sel="(all) and (name H*)",
            acceptors_sel="(all) and (name N* O*)",
        )
        # the reloaded reader can order bonds differently
        assert set(map(tuple, edges(obj))) == set(map(tuple, expected))
    assert mol.calculations._hbonds._hba_trajectory is mol.universe.trajectory


def test_hbonds_remove(fetch):
    mol = fetch("2M6Q")
    obj = mol.calculations.hbonds().object
    n = len(obj.data.edges)
    del mol.calculations["hbonds"]
    assert "hbonds" not in mol.frame_attributes
    assert not mol.calculations._hbonds.is_setup
    # the data object keeps the bonds it last had
    mol.set_frame(3)
    assert len(obj.data.edges) == n


def test_hbonds_pickle(universe):
    mol = mn.Molecule(universe, name="TestPickleHBonds")
    mol.calculations.hbonds("resid 1-50", updating=False, d_a_cutoff=3.2)
    restored = pickle.loads(pickle.dumps(mol))
    # the function isn't pickled, the saved settings are
    assert "hbonds" not in restored.calculations
    restored.set_frame(1)
    hbonds = restored.calculations._hbonds
    assert hbonds.is_setup
    assert hbonds.cutoffs["d_a_cutoff"] == pytest.approx(3.2)
    assert restored.props.hbonds.selection == "resid 1-50"
    assert not restored.calculations.updating("hbonds")
    assert hbonds.object == mol.calculations._hbonds.object


def test_style_hbonds_follows_atoms(fetch):
    mol = fetch("2M6Q")
    obj = mol.calculations.hbonds().object
    pairs = edges(obj)
    offset = np.array([1.0, 2.0, 3.0])
    with mol.tree.reset() as (atoms, join):
        (
            atoms
            >> g.SetPosition(offset=tuple(offset))
            >> StyleHBonds(data_object=obj, dashes=1, quality=2)
            >> join
        )
    position = mol.named_attribute("position", evaluate=True)
    # quality 2 gives cylinders with two rings of 8 vertices, one per hydrogen bond
    cylinders = position.reshape(len(pairs), 16, 3)
    atom_position = mol.position + offset
    midpoints = atom_position[pairs].mean(axis=1)
    np.testing.assert_allclose(cylinders.mean(axis=1), midpoints, atol=1e-5)
    # dashes are evenly split along the bond
    with mol.tree.reset() as (atoms, join):
        atoms >> StyleHBonds(data_object=obj, dashes=3, quality=2) >> join
    assert evaluated_vertex_count(mol.object) == len(pairs) * 3 * 16


def test_style_hbonds_selection(fetch):
    mol = fetch("2M6Q")
    obj = mol.calculations.hbonds().object
    pairs = edges(obj)
    is_backbone = mol.named_attribute("is_backbone")
    n_expected = np.all(is_backbone[pairs], axis=1).sum()
    assert 0 < n_expected < len(pairs)
    with mol.tree.reset() as (atoms, join):
        (
            atoms
            >> StyleHBonds(
                selection=g.NamedAttribute.boolean("is_backbone"),
                data_object=obj,
                dashes=1,
                quality=2,
            )
            >> join
        )
    assert evaluated_vertex_count(mol.object) == n_expected * 16


@pytest.mark.parametrize("distance, angle", [(3.0, 150.0), (3.5, 130.0)])
def test_find_hbonds_node(fetch, distance, angle):
    mol = fetch("2M6Q")
    with mol.tree.reset() as (atoms, join):
        find = FindHBonds(atoms=atoms, distance=distance, angle=np.radians(angle))
        (
            atoms
            >> g.StoreNamedAttribute.point.integer(
                name="acceptor", value=find.o.acceptor_index
            )
            >> g.StoreNamedAttribute.point.boolean(
                name="is_hbonded", value=find.o.is_hbonded
            )
            >> join
        )
    acceptor = mol.named_attribute("acceptor", evaluate=True)
    is_hbonded = mol.named_attribute("is_hbonded", evaluate=True)
    np.testing.assert_array_equal(is_hbonded, acceptor >= 0)
    found = {(h, a) for h, a in enumerate(acceptor) if a >= 0}
    expected = {
        tuple(pair)
        for pair in expected_pairs(
            mol.universe.copy(),
            0,
            d_a_cutoff=distance,
            d_h_a_angle_cutoff=angle,
            donors_sel="element N O",
            hydrogens_sel="element H D",
            acceptors_sel="element N O",
        )
    }
    assert len(expected) > 0
    # each hydrogen bonds to at most one acceptor in nodes, so a hydrogen
    # bonding to two (bifurcated) gives one of them
    assert found <= expected
    bifurcated = len(expected) - len({h for h, _ in expected})
    assert len(found) >= len(expected) - bifurcated


def test_style_hbonds_find_hbonds_source(fetch):
    mol = fetch("2M6Q")
    obj = mol.calculations.hbonds().object
    pairs = edges(obj)
    offset = np.array([1.0, 2.0, 3.0])
    with mol.tree.reset() as (atoms, join):
        (
            atoms
            >> g.SetPosition(offset=tuple(offset))
            >> StyleHBonds(source="Find HBonds", dashes=1, quality=2)
            >> join
        )
    position = mol.named_attribute("position", evaluate=True)
    # the same bonds as MDAnalysis at the default cutoffs, following the atoms
    cylinders = position.reshape(len(pairs), 16, 3)
    midpoints = (mol.position + offset)[pairs].mean(axis=1)
    np.testing.assert_allclose(
        np.sort(cylinders.mean(axis=1), axis=0), np.sort(midpoints, axis=0), atol=1e-5
    )
    # a bond is shown only when both atoms are selected
    is_backbone = mol.named_attribute("is_backbone")
    n_expected = np.all(is_backbone[pairs], axis=1).sum()
    with mol.tree.reset() as (atoms, join):
        (
            atoms
            >> StyleHBonds(
                source="Find HBonds",
                selection=g.NamedAttribute.boolean("is_backbone"),
                dashes=1,
                quality=2,
            )
            >> join
        )
    assert evaluated_vertex_count(mol.object) == n_expected * 16


def test_hbonds_init_operator(fetch):
    mol = fetch("2M6Q")
    hbonds = mol.calculations._hbonds
    assert not hbonds.is_setup
    mn_ops.hbonds_init(uuid=mol.uuid, selection="protein")
    assert hbonds.is_setup
    assert mol.props.hbonds.selection == "protein"
    n = len(hbonds.object.data.edges)
    # a second run updates the selection
    mn_ops.hbonds_init(uuid=mol.uuid, selection="resid 1-20")
    assert mol.props.hbonds.selection == "resid 1-20"
    assert len(hbonds.object.data.edges) < n
