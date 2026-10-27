import bpy
import numpy as np
import pytest
from mathutils import Euler, Matrix
import molecularnodes as mn
from molecularnodes.ui.ops import MN_OT_Align
from .constants import data_dir
from .utils import mn_ops


def _random_matrix(seed: int) -> Matrix:
    rng = np.random.default_rng(seed)
    rotation = Euler(rng.uniform(0, 2 * np.pi, 3).tolist()).to_matrix().to_4x4()
    return Matrix.Translation(rng.uniform(-5, 5, 3).tolist()) @ rotation


def _world_positions(mol: mn.Molecule, selection: str) -> np.ndarray:
    atoms = mol.universe.select_atoms(selection)
    local = atoms.positions * mol.world_scale
    matrix = np.array(mol.object.matrix_world)
    return local @ matrix[:3, :3].T + matrix[:3, 3]


def _homolog_pair() -> tuple[mn.Molecule, mn.Molecule]:
    """Two copies of the same structure with differently transformed objects;
    the mobile copy has point mutations, so the sequences differ."""
    fixed = mn.Molecule.load(data_dir / "1l58.cif")
    mobile = mn.Molecule.load(data_dir / "1l58.cif")
    mobile.atoms.residues[20:30].resnames = "GLY"
    fixed.object.matrix_world = _random_matrix(1)
    mobile.object.matrix_world = _random_matrix(2)
    return fixed, mobile


@pytest.mark.parametrize("method", ["sequence", "structure"])
def test_align_recovers_transform(method):
    fixed, mobile = _homolog_pair()

    # a missing loop on the mobile side has to be bridged by the alignment
    result = mobile.align_to(
        fixed, method=method, selection="protein and not resid 50:60"
    )

    assert result.rmsd < 0.01
    # outlier removal is quantile based, so on an exact copy (all distances ~0)
    # it prunes many anchors without changing the fit
    assert result.n_anchors >= 20
    ca = "protein and name CA"
    assert np.allclose(
        _world_positions(mobile, ca), _world_positions(fixed, ca), atol=1e-3
    )


def test_align_select_chains():
    fixed = mn.Molecule.load(data_dir / "8H1B.cif")
    mobile = mn.Molecule.load(data_dir / "8H1B.cif")
    mobile.object.matrix_world = _random_matrix(3)
    result = mobile.align_to(fixed, selection="chainID B", fixed_selection="chainID B")
    assert not result.merged_chains
    assert result.rmsd < 0.01
    ca = "protein and name CA"
    assert np.allclose(
        _world_positions(mobile, ca), _world_positions(fixed, ca), atol=1e-3
    )


@pytest.mark.parametrize("method", ["sequence", "structure"])
def test_align_chain_count_mismatch(method):
    # like a Fab's separate heavy and light chains against a single-chain scFv:
    # the mobile copy is split into two chains, the fixed one is a single chain
    fixed, mobile = _homolog_pair()
    mobile.atoms[mobile.atoms.resids > 80].chainIDs = "B"

    result = mobile.align_to(fixed, method=method)

    assert result.merged_chains
    assert result.rmsd < 0.01
    ca = "protein and name CA"
    assert np.allclose(
        _world_positions(mobile, ca), _world_positions(fixed, ca), atol=1e-3
    )


def test_align_onto_itself():
    mol = mn.Molecule.load(data_dir / "1l58.cif")
    with pytest.raises(ValueError, match="itself"):
        mol.align_to(mol)


def test_op_align():
    fixed, mobile = _homolog_pair()
    bpy.ops.object.select_all(action="DESELECT")
    mobile.object.select_set(True)
    fixed.object.select_set(True)
    view_layer = bpy.context.view_layer
    assert view_layer is not None
    view_layer.objects.active = fixed.object

    assert mn_ops.align() == {"FINISHED"}
    ca = "protein and name CA"
    assert np.allclose(
        _world_positions(mobile, ca), _world_positions(fixed, ca), atol=1e-3
    )


def test_op_align_chain_summary():
    # the dialog lists each molecule's chains to help write the selections
    mol = mn.Molecule.load(data_dir / "8H1B.cif")
    summary = MN_OT_Align._chain_summary(mol)
    assert [part.split()[0] for part in summary.split(", ")] == ["A", "B", "C", "D"]
