"""
Superimpose one molecule onto another.

The fit is done with biotite on the atoms in Å, in each molecule's own local
frame, and applied to the mobile molecule's Blender object as its world matrix.
The coordinates themselves are left untouched, so the alignment can be undone,
moved by hand afterwards and survives trajectory frame changes and reloads.
"""

from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal
import biotite.structure as struc
import MDAnalysis as mda
import numpy as np
from mathutils import Matrix

if TYPE_CHECKING:
    from .base import Molecule

ALIGN_METHODS = ("sequence", "structure")

# the atoms considered by default: superimpose_homologs anchors on CA / P atoms,
# superimpose_structural_homologs only works on peptide chains
DEFAULT_SELECTIONS = {"sequence": "protein or nucleic", "structure": "protein"}


@dataclass
class AlignResult:
    """The outcome of :func:`align`.

    Attributes
    ----------
    matrix : np.ndarray
        The 4x4 transform, in Å, that maps the mobile atoms onto the fixed atoms,
        both in their objects' local frames.
    rmsd : float
        RMSD in Å between the anchor atoms after the fit.
    n_anchors : int
        Number of atom pairs the final fit was made on.
    """

    matrix: np.ndarray
    rmsd: float
    n_anchors: int


def atom_array_from_atoms(atoms: mda.AtomGroup) -> struc.AtomArray:
    """Build a biotite ``AtomArray`` from an MDAnalysis ``AtomGroup``, at the
    universe's current frame.

    Only the annotations biotite's superimposition needs are filled. Chains come
    from ``chainIDs``, falling back to ``segids`` for topologies without them.
    """
    array = struc.AtomArray(len(atoms))
    array.coord = atoms.positions.astype(np.float32)
    array.res_id = atoms.resids
    array.res_name = atoms.resnames
    array.atom_name = atoms.names
    if hasattr(atoms, "elements"):
        array.element = atoms.elements
    if hasattr(atoms, "chainIDs"):
        array.chain_id = atoms.chainIDs
    else:
        array.chain_id = atoms.segids
    if hasattr(atoms, "icodes"):
        array.ins_code = atoms.icodes
    return array


def _select(mol: "Molecule", selection: str, role: str) -> struc.AtomArray:
    atoms = mol.universe.select_atoms(selection)
    if len(atoms) == 0:
        raise ValueError(f"The {role} selection '{selection}' matches no atoms")
    return atom_array_from_atoms(atoms)


def align(
    mobile: "Molecule",
    fixed: "Molecule",
    method: Literal["sequence", "structure"] = "sequence",
    selection: str | None = None,
    fixed_selection: str | None = None,
) -> AlignResult:
    """Superimpose ``mobile`` onto ``fixed`` by setting ``mobile``'s object
    transform.

    Parameters
    ----------
    mobile : Molecule
        The molecule to move.
    fixed : Molecule
        The reference molecule, which stays where it is.
    method : {"sequence", "structure"}, optional
        ``"sequence"`` (default) pairs residues by sequence alignment and fits on
        their CA or P atoms, iteratively dropping outliers, with
        :func:`biotite.structure.superimpose_homologs`. It copes with mutations,
        insertions and missing residues. ``"structure"`` pairs residues by a
        structural alphabet (3Di) with
        :func:`biotite.structure.superimpose_structural_homologs`, for remote
        homologs with low sequence identity; proteins only.
    selection : str | None, optional
        MDAnalysis selection of the mobile atoms to fit on. Defaults to
        ``"protein or nucleic"`` for ``"sequence"`` and ``"protein"`` for
        ``"structure"``.
    fixed_selection : str | None, optional
        MDAnalysis selection of the fixed atoms to fit on, with the same default.

    The two selections must contain the same number of chains, in matching
    order; select single chains (e.g. ``"chainID A"`` and ``"chainID H"``) to
    align chains that are labelled differently or to align a subunit.

    Returns
    -------
    AlignResult
        The fitted transform, RMSD and number of anchor atoms.
    """
    if method not in ALIGN_METHODS:
        raise ValueError(f"method must be one of {ALIGN_METHODS}, not '{method}'")
    if mobile.object == fixed.object:
        raise ValueError("Cannot align a molecule onto itself")

    default = DEFAULT_SELECTIONS[method]
    fixed_array = _select(fixed, fixed_selection or default, "fixed")
    mobile_array = _select(mobile, selection or default, "mobile")

    n_fixed = struc.get_chain_count(fixed_array)
    n_mobile = struc.get_chain_count(mobile_array)
    if n_fixed != n_mobile:
        raise ValueError(
            f"The fixed selection has {n_fixed} chains and the mobile selection"
            f" {n_mobile}; select matching chains on each side"
            " (e.g. 'chainID A' and 'chainID H')"
        )

    if method == "sequence":
        fitted, transform, fixed_idx, mobile_idx = struc.superimpose_homologs(
            fixed_array, mobile_array
        )
    else:
        fitted, transform, fixed_idx, mobile_idx = (
            struc.superimpose_structural_homologs(fixed_array, mobile_array)
        )

    if len(fixed_idx) == 0:
        raise ValueError("No corresponding atoms were found to align on")
    rmsd = float(struc.rmsd(fixed_array.coord[fixed_idx], fitted.coord[mobile_idx]))
    matrix = transform.as_matrix()
    if matrix.ndim == 3:
        matrix = matrix[0]

    # the fit is in Å in each object's local frame; convert it to map the mobile
    # object's local units onto the fixed object's, then place the mobile object
    # relative to the fixed one
    local = np.array(matrix, dtype=float)
    local[:3, :3] *= fixed.world_scale / mobile.world_scale
    local[:3, 3] *= fixed.world_scale
    mobile.object.matrix_world = fixed.object.matrix_world @ Matrix(local.tolist())

    return AlignResult(matrix=matrix, rmsd=rmsd, n_anchors=len(fixed_idx))
