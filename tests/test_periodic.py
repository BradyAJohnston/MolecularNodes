import itertools
from typing import cast
import bpy
import MDAnalysis as mda
import numpy as np
import pytest
from MDAnalysis.lib.distances import apply_PBC, minimize_vectors
from MDAnalysis.lib.mdamath import triclinic_vectors
from MDAnalysis.tests.datafiles import DCD, GRO, PSF, XTC
from nodebpy.nodes.geometry import CurveToMesh, SetPosition, StoreNamedAttribute
from scipy.spatial import Voronoi
import molecularnodes as mn
from molecularnodes.nodes.geometry import (
    FragmentID,
    PeriodicCell,
    PeriodicOffset,
    PeriodicUnwrap,
)
from .constants import data_dir

BOXES = {
    "cubic": [85.5, 85.5, 85.5, 90.0, 90.0, 90.0],
    "orthorhombic": [53.8, 44.0, 57.2, 90.0, 90.0, 90.0],
    "hexagonal": [98.0, 98.0, 130.7, 90.0, 90.0, 120.0],
    "dodecahedron": [82.68, 82.68, 82.68, 60.0, 60.0, 90.0],
    "dodecahedron-hexagon": [80.0, 80.0, 80.0, 60.0, 60.0, 60.0],
    "truncated-octahedron": [80.0, 80.0, 80.0, 70.53, 109.47, 70.53],
}
BOX_NAMES = ["a", "b", "c", "alpha", "beta", "gamma"]
# the default world scale of a Molecule, 1 Å is 0.1 world units
WORLD_SCALE = 0.1
# Å, and the tolerance for comparing positions read back from the mesh
ATOL = 1e-2
# fractional offset that splits molecules across the box before unwrapping
SHIFT = np.array([0.37, 0.21, 0.5])


def _world_box(box) -> np.ndarray:
    """Box dimensions with the lengths in world units, as the nodes take them"""
    box = np.asarray(box, dtype=float)
    return np.concatenate([box[:3] * WORLD_SCALE, box[3:]])


def _box_inputs(box) -> dict:
    return {"update": False, **dict(zip(BOX_NAMES, map(float, _world_box(box))))}


def _box_values(node) -> list[float]:
    """The box lengths and angles set on a node's inputs"""
    return [
        cast(bpy.types.NodeSocketFloat, node.node.inputs[name]).default_value
        for name in BOX_NAMES
    ]


def _evaluated_positions(mol: mn.Molecule) -> np.ndarray:
    """Positions after the node tree, in Å"""
    return mol.named_attribute("position", evaluate=True) / mol.world_scale


def _assert_inside_box(positions: np.ndarray, box) -> None:
    fractional = positions @ np.linalg.inv(triclinic_vectors(box))
    assert fractional.min() > -1e-4 and fractional.max() < 1 + 1e-4


def _assert_same_image(result: np.ndarray, expected: np.ndarray, box) -> None:
    """Equal, allowing for atoms on a box face wrapped to the opposite one"""
    fractional = (result - expected) @ np.linalg.inv(triclinic_vectors(box))
    fractional -= np.round(fractional)
    assert np.abs(fractional @ triclinic_vectors(box)).max() < ATOL


def _assert_same_up_to_lattice(result, expected, groups, box) -> None:
    """Each group is the same shape as expected, moved by a whole box vector"""
    _, first, inverse = np.unique(groups, return_index=True, return_inverse=True)
    diff = result - expected
    shape_error = diff - diff[first][inverse]
    assert np.abs(shape_error).max() < ATOL
    lattice = diff[first] @ np.linalg.inv(triclinic_vectors(box))
    assert np.allclose(lattice, np.round(lattice), atol=1e-3)


@pytest.fixture
def mol_triclinic():
    # rhombic dodecahedron with bonds, whole molecules and 11779 fragments
    u = mda.Universe(data_dir / "md_ppr/md.tpr", data_dir / "md_ppr/md.gro")
    return mn.Molecule(u)


@pytest.mark.parametrize("box_name", BOXES)
@pytest.mark.parametrize(
    "offset_type, offset",
    [
        ("Factor", (0.0, 0.0, 0.0)),
        ("Factor", (0.37, -0.21, 1.5)),
        ("Distance", (1.3, -2.2, 0.4)),
    ],
)
def test_periodic_offset(mol_triclinic, box_name, offset_type, offset):
    mol = mol_triclinic
    box = np.array(BOXES[box_name])
    positions = mol.named_attribute("position") / mol.world_scale
    with mol.tree.reset() as (atoms, join):
        node = PeriodicOffset(
            offset_type=offset_type,
            offset_factor=offset,
            offset_distance=offset,
            **_box_inputs(box),
        )
        atoms >> SetPosition(position=node) >> join
    if offset_type == "Factor":
        shift = np.array(offset) @ triclinic_vectors(box)
    else:
        # world units, where 1 unit is 1 nm
        shift = np.array(offset) / mol.world_scale
    expected = apply_PBC((positions + shift).astype(np.float32), box)
    result = _evaluated_positions(mol)
    _assert_inside_box(result, box)
    _assert_same_image(result, expected, box)


@pytest.mark.parametrize("box_name", BOXES)
def test_periodic_offset_compact(mol_triclinic, box_name):
    # each position at its image nearest the box centre, as GROMACS -ur compact
    mol = mol_triclinic
    box = np.array(BOXES[box_name])
    offset = np.array([0.37, -0.21, 1.5])
    positions = mol.named_attribute("position") / mol.world_scale
    with mol.tree.reset() as (atoms, join):
        node = PeriodicOffset(
            offset_factor=tuple(offset), wrap="Compact", **_box_inputs(box)
        )
        atoms >> SetPosition(position=node) >> join
    vectors = triclinic_vectors(box)
    centre = vectors.sum(axis=0) / 2
    shifted = positions + offset @ vectors
    expected = centre + minimize_vectors((shifted - centre).astype(np.float32), box)
    result = _evaluated_positions(mol)
    _assert_same_image(result, expected, box)
    # no farther from the centre than the nearest image, allowing for ties
    distance = np.linalg.norm(result - centre, axis=1)
    assert np.all(distance <= np.linalg.norm(expected - centre, axis=1) + ATOL)


@pytest.mark.parametrize("node", [PeriodicOffset, PeriodicUnwrap])
def test_periodic_no_box_unchanged(mol_triclinic, node):
    mol = mol_triclinic
    positions = mol.named_attribute("position")
    with mol.tree.reset() as (atoms, join):
        atoms >> SetPosition(position=node(update=False)) >> join
    np.testing.assert_allclose(
        mol.named_attribute("position", evaluate=True), positions, atol=1e-6
    )


def test_fragment_id(mol_triclinic):
    mol = mol_triclinic
    expected = mol.universe.atoms.fragindices
    np.testing.assert_array_equal(mol.named_attribute("fragment_id"), expected)
    # read by the Fragment ID node
    with mol.tree.reset() as (atoms, join):
        store = StoreNamedAttribute.point.integer(
            name="node_fragment_id", value=FragmentID()
        )
        atoms >> store >> join
    np.testing.assert_array_equal(
        mol.named_attribute("node_fragment_id", evaluate=True), expected
    )
    # a topology without bonds has no fragments
    assert "fragment_id" not in mn.Molecule(mda.Universe(GRO)).list_attributes()


def _break_and_unwrap(mol, box, place_groups, group_id=lambda: None):
    """Split molecules across the box with Periodic Offset, then unwrap them"""
    with mol.tree.reset() as (atoms, join):
        offset = PeriodicOffset(offset_factor=tuple(SHIFT), **_box_inputs(box))
        unwrap = PeriodicUnwrap(
            group_id=group_id(), place_groups=place_groups, **_box_inputs(box)
        )
        (atoms >> SetPosition(position=offset) >> SetPosition(position=unwrap) >> join)
    positions = mol.named_attribute("position") / mol.world_scale
    broken = apply_PBC(
        (positions + SHIFT @ triclinic_vectors(box)).astype(np.float32),
        box,
    )
    return broken, _evaluated_positions(mol)


@pytest.mark.parametrize("place_groups", [False, True])
def test_periodic_unwrap(mol_triclinic, place_groups):
    mol = mol_triclinic
    u = mol.universe
    box = u.dimensions
    positions = mol.named_attribute("position") / mol.world_scale
    broken, result = _break_and_unwrap(mol, box, place_groups, FragmentID)
    # the offset splits molecules, which the unwrap puts back together
    assert u.bonds.values().max() < 2.0
    u.atoms.positions = broken
    assert u.bonds.values().max() > 10.0
    u.atoms.positions = result
    assert u.bonds.values().max() < 2.0
    fragments = u.atoms.fragindices
    shifted = positions + SHIFT @ triclinic_vectors(box)
    _assert_same_up_to_lattice(result, shifted, fragments, box)
    if place_groups:
        u.atoms.positions = result
        centres = u.atoms.center_of_geometry(compound="fragments")
        fractional = centres @ np.linalg.inv(triclinic_vectors(box))
        assert fractional.min() > -0.1 and fractional.max() < 1.1


def _unwrap(mol, box, **kwargs):
    """Unwrap by fragment, returning positions in Å"""
    with mol.tree.reset() as (atoms, join):
        node = PeriodicUnwrap(group_id=FragmentID(), **kwargs, **_box_inputs(box))
        atoms >> SetPosition(position=node) >> join
    return _evaluated_positions(mol)


@pytest.mark.parametrize("box_name", BOXES)
def test_periodic_unwrap_whole_stays_in_place(mol_triclinic, box_name):
    # molecules that are already whole are left where they are without Place
    # Groups, in any box, as each centre moves to its image nearest the average
    mol = mol_triclinic
    box = np.array(BOXES[box_name])
    positions = mol.named_attribute("position") / mol.world_scale
    result = _unwrap(mol, box, place_groups=False)
    # only groups shorter than half the box can be unwrapped, which excludes
    # the protein in the orthorhombic box, smaller than the protein itself
    half = np.diag(triclinic_vectors(box)).min() / 2
    fits = np.concatenate(
        [
            np.full(len(f), np.ptp(f.positions, axis=0).max() < half)
            for f in mol.universe.atoms.fragments
        ]
    )[np.argsort(np.concatenate([f.indices for f in mol.universe.atoms.fragments]))]
    assert fits.mean() > 0.8
    assert np.abs(result - positions)[fits].max() < ATOL


@pytest.mark.parametrize("box_name", BOXES)
def test_periodic_unwrap_place_groups_compact(mol_triclinic, box_name):
    # each group centre is at its image nearest the box centre
    mol = mol_triclinic
    box = np.array(BOXES[box_name])
    result = _unwrap(mol, box, place_groups=True, wrap="Compact")
    vectors = triclinic_vectors(box)
    groups = mol.universe.atoms.fragindices
    count = np.bincount(groups)[:, None]
    # the circular mean centre the node places, recovered from the result: the
    # image of it nearest each group's plain average
    angle = 2 * np.pi * (result @ np.linalg.inv(vectors))
    mean_cos = np.stack([np.bincount(groups, np.cos(angle[:, i])) for i in range(3)], 1)
    mean_sin = np.stack([np.bincount(groups, np.sin(angle[:, i])) for i in range(3)], 1)
    circular = (np.arctan2(mean_sin, mean_cos) / (2 * np.pi)) @ vectors
    average = np.stack([np.bincount(groups, result[:, i]) for i in range(3)], 1) / count
    placed = average + minimize_vectors((circular - average).astype(np.float32), box)
    box_centre = vectors.sum(axis=0) / 2
    nearest = box_centre + minimize_vectors(
        (placed - box_centre).astype(np.float32), box
    )
    # equal, or an equally near image on a cell face
    distance = np.linalg.norm(placed - box_centre, axis=1)
    assert np.all(distance <= np.linalg.norm(nearest - box_centre, axis=1) + ATOL)
    _assert_same_image(placed, nearest, box)


@pytest.mark.parametrize("place_groups", [False, True])
def test_periodic_unwrap_default_group_per_atom(mol_triclinic, place_groups):
    # unlinked, Group ID is the ID or index, so each atom is its own group
    mol = mol_triclinic
    box = mol.universe.dimensions
    broken, result = _break_and_unwrap(mol, box, place_groups)
    _assert_same_image(result, broken, box)
    if place_groups:
        _assert_inside_box(result, box)
    else:
        assert np.abs(result - broken).max() < ATOL


def test_periodic_box_updates():
    mol = mn.Molecule.load(GRO, XTC)
    with mol.tree.reset() as (atoms, join):
        offset = PeriodicOffset()
        unwrap = PeriodicUnwrap()
        atoms >> SetPosition(position=offset) >> SetPosition(position=unwrap) >> join
        cell = PeriodicCell()
        cell.o.mesh >> join
    for frame in [1, 5]:
        mol.set_frame(frame)
        dimensions = mol.universe.trajectory.ts.dimensions
        for node in [offset, unwrap, cell]:
            np.testing.assert_allclose(
                _box_values(node), _world_box(dimensions), rtol=1e-6
            )


def test_periodic_box_updates_without_dimensions():
    mol = mn.Molecule.load(PSF, DCD)
    with mol.tree.reset() as (atoms, join):
        node = PeriodicOffset()
        atoms >> SetPosition(position=node) >> join
    mol.set_frame(1)
    assert _box_values(node) == [0] * 6


def _dataset_universe(name: str) -> mda.Universe:
    """A universe from MDAnalysisData, skipping when it isn't downloaded"""
    from MDAnalysisData import datasets

    try:
        data = getattr(datasets, f"fetch_{name}")(download_if_missing=False)
    except OSError:
        pytest.skip(f"MDAnalysisData '{name}' is not downloaded")
    return mda.Universe(data.topology, data.trajectory)


@pytest.mark.parametrize(
    "name",
    [
        "adk_equilibrium",
        "PEG_1chain",
        "membrane_peptide",
        "yiip_equilibrium_short",
        "ifabp_water",
    ],
)
@pytest.mark.parametrize("place_groups", [False, True])
def test_periodic_unwrap_datasets(name, place_groups):
    u = _dataset_universe(name)
    if name == "yiip_equilibrium_short":
        u.guess_TopologyAttrs(to_guess=["bonds"], vdwradii={"Z": 1.39, "D": 1.2})
    mol = mn.Molecule(u)
    box = u.dimensions
    if box is None:
        # no box, as in the fitted ifabp_water trajectory: nothing changes
        positions = mol.named_attribute("position")
        with mol.tree.reset() as (atoms, join):
            node = PeriodicUnwrap(group_id=FragmentID())
            atoms >> SetPosition(position=node) >> join
        np.testing.assert_allclose(
            mol.named_attribute("position", evaluate=True), positions, atol=1e-6
        )
        return
    broken, result = _break_and_unwrap(mol, box, place_groups, FragmentID)
    atoms = u.atoms
    assert atoms is not None
    atoms.positions = broken
    expected = atoms.unwrap(compound="fragments")
    _assert_same_up_to_lattice(result, expected, atoms.fragindices, box)


def _cell(mol, box, wrap, output):
    """Vertices (Å) and edge count of a Periodic Cell output, as a mesh"""
    with mol.tree.reset() as (atoms, join):
        node = PeriodicCell(wrap=wrap, **_box_inputs(box))
        if output == "edges":
            node.o.edges >> CurveToMesh() >> join
        else:
            node.o.mesh >> join
    obj = mol.object.evaluated_get(bpy.context.evaluated_depsgraph_get())
    mesh = cast(bpy.types.Mesh, obj.data)
    vertices = np.array([v.co for v in mesh.vertices]) / mol.world_scale
    return vertices, len(mesh.edges)


def _same_points(a, b, atol=ATOL):
    assert len(a) == len(b)
    distance = np.linalg.norm(a[:, None] - b[None], axis=-1)
    assert distance.min(axis=1).max() < atol


# corners and edges of each compact cell: a box, a hexagonal prism, a rhombic
# dodecahedron or a truncated octahedron
COMPACT_CELLS = {
    "cubic": (8, 12),
    "orthorhombic": (8, 12),
    "hexagonal": (12, 18),
    "dodecahedron": (14, 24),
    "dodecahedron-hexagon": (14, 24),
    "truncated-octahedron": (24, 36),
}


@pytest.mark.parametrize("box_name", BOXES)
def test_periodic_cell_box(mol_triclinic, box_name):
    box = np.array(BOXES[box_name])
    vectors = triclinic_vectors(box)
    corners = np.array(list(itertools.product((0, 1), repeat=3))) @ vectors
    vertices, _ = _cell(mol_triclinic, box, "Box", "mesh")
    _same_points(vertices, corners)
    _, n_edges = _cell(mol_triclinic, box, "Box", "edges")
    assert n_edges == 12


@pytest.mark.parametrize("box_name", BOXES)
def test_periodic_cell_compact(mol_triclinic, box_name):
    # the Wigner-Seitz cell around the box centre, as the simulation box
    # annotation draws it with scipy
    box = np.array(BOXES[box_name])
    vectors = triclinic_vectors(box).astype(np.float64)
    shifts = np.array(list(itertools.product((0, -1, 1), repeat=3))) @ vectors
    voronoi = Voronoi(shifts)
    region = voronoi.regions[voronoi.point_region[0]]
    expected = voronoi.vertices[region] + vectors.sum(axis=0) / 2
    # where more than four cells meet, as at the corners of a rhombic
    # dodecahedron, qhull gives several vertices at one corner
    expected = expected[
        [
            i
            for i in range(len(expected))
            if np.all(np.linalg.norm(expected[:i] - expected[i], axis=1) > ATOL)
        ]
    ]
    n_corners, n_edges = COMPACT_CELLS[box_name]
    vertices, _ = _cell(mol_triclinic, box, "Compact", "mesh")
    assert len(vertices) == n_corners
    _same_points(vertices, expected)
    _, edges = _cell(mol_triclinic, box, "Compact", "edges")
    assert edges == n_edges


def test_periodic_cell_no_box(mol_triclinic):
    vertices, n_edges = _cell(mol_triclinic, np.zeros(6), "Compact", "mesh")
    assert len(vertices) == 0 and n_edges == 0
