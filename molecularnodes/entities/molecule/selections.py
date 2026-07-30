import numpy as np
from biotite.structure import AtomArray, AtomArrayStack, filter


def select_amino_acids(arr: AtomArray | AtomArrayStack) -> np.ndarray:
    return filter.filter_amino_acids(arr)


def select_atom_names(
    arr: AtomArray | AtomArrayStack,
    atom_name: str | list[str] | tuple[str, ...] | np.ndarray,
) -> np.ndarray:
    return np.isin(arr.get_annotation("atom_name"), atom_name)


def select_canonical_amino_acids(arr: AtomArray | AtomArrayStack) -> np.ndarray:
    return filter.filter_canonical_amino_acids(arr)


def select_canonical_nucleotides(arr: AtomArray | AtomArrayStack) -> np.ndarray:
    return filter.filter_canonical_nucleotides(arr)
