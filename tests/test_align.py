import bpy
import numpy as np
import pytest
from mathutils import Euler, Matrix
import molecularnodes as mn
from .constants import data_dir


def _random_matrix(seed: int) -> Matrix:
    rng = np.random.default_rng(seed)
    rotation = Euler(rng.uniform(0, 2 * np.pi, 3)).to_matrix().to_4x4()
    return Matrix.Translation(rng.uniform(-5, 5, 3)) @ rotation


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
    mobile.universe.residues[20:30].resnames = "GLY"
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


def test_align_chain_count_mismatch():
    fixed = mn.Molecule.load(data_dir / "8H1B.cif")
    mobile = mn.Molecule.load(data_dir / "8H1B.cif")
    with pytest.raises(ValueError, match="chains"):
        mobile.align_to(fixed, selection="chainID A")

    # selecting a single chain on each side makes them comparable
    result = mobile.align_to(fixed, selection="chainID A", fixed_selection="chainID A")
    assert result.rmsd < 0.01


def test_align_onto_itself():
    mol = mn.Molecule.load(data_dir / "1l58.cif")
    with pytest.raises(ValueError, match="itself"):
        mol.align_to(mol)


def test_op_align():
    fixed, mobile = _homolog_pair()
    bpy.ops.object.select_all(action="DESELECT")
    mobile.object.select_set(True)
    fixed.object.select_set(True)
    bpy.context.view_layer.objects.active = fixed.object

    assert bpy.ops.mn.align() == {"FINISHED"}
    ca = "protein and name CA"
    assert np.allclose(
        _world_positions(mobile, ca), _world_positions(fixed, ca), atol=1e-3
    )
