import numpy as np
from biotite.structure import filter, AtomArray
from typing import Callable


class StructureSelector:
    def __init__(self):
        self.mask: np.ndarray | None = None
        self.pending_selections: list[Callable] = []

    def _update_mask(self, operation):
        self.pending_selections.append(operation)
        return self

    def clear_selections(self):
        self.pending_selections = []

    def evaluate_on_array(self, array: AtomArray) -> np.ndarray:
        """Evaluates all pending operations on the given array"""
        mask = np.ones(array.array_length(), dtype=bool)
        for operation in self.pending_selections:
            mask = np.logical_and(mask, operation(array))
        self.mask = mask
        return mask

    def get_selection(self, array):
        """Returns the structure filtered by the current mask"""
        mask = self.evaluate_on_array(array)
        if mask is None:
            return array
        return array[mask]

    def amino_acids(self):
        return self._update_mask(select_amino_acids)

    def atom_name(self, atom_name: list[str] | tuple[str, ...] | np.ndarray):
        return self._update_mask(lambda arr: select_atom_names(arr, atom_name))

    def is_canonical_amino_acid(self):
        return self._update_mask(select_canonical_amino_acids)

    def is_canonical_nucleotide(self):
        return self._update_mask(select_canonical_nucleotides)

    def is_carbohydrate(self):
        return self._update_mask(select_carbohydrates)

    def chain_id(self, chain_id: list[str] | tuple[str, ...] | np.ndarray):
        return self._update_mask(lambda arr: select_chain_id(arr, chain_id))

    def element(self, element: list[str] | tuple[str, ...] | np.ndarray):
        return self._update_mask(lambda arr: select_element(arr, element))

    # def first_altloc(self):
    #     return self._update_mask(select_first_altloc)

    def is_hetero(self):
        return self._update_mask(select_hetero)

    # def highest_occupancy_altloc(self):
    #     return self._update_mask(select_highest_occupancy_altloc)

    # def inscode(self, inscode):
    #     return self._update_mask(lambda arr: select_inscode(arr, inscode))

    # def intersection(self):
    #     return self._update_mask(select_intersection)

    def is_ligand(self):
        return self._update_mask(select_ligand)

    def linear_bond_continuity(self):
        return self._update_mask(select_linear_bond_continuity)

    def is_monoatomic_ion(self):
        return self._update_mask(select_monoatomic_ions)

    def is_nucleotide(self):
        return self._update_mask(select_nucleotides)

    def peptide_backbone(self):
        return self._update_mask(select_peptide_backbone)

    def phosphate_backbone(self):
        return self._update_mask(select_phosphate_backbone)

    def polymer(self):
        return self._update_mask(select_polymer)

    def res_id(self, num):
        return self._update_mask(lambda arr: select_res_id(arr, num))

    def res_ids(self, nums):
        return self._update_mask(lambda arr: select_res_ids(arr, nums))

    def res_name(self, res_name):
        return self._update_mask(lambda arr: select_res_name(arr, res_name))

    def solvent(self):
        return self._update_mask(select_solvent)

    def not_amino_acids(self):
        return self._update_mask(lambda arr: ~select_amino_acids(arr))

    def not_atom_names(self, atomname):
        return self._update_mask(lambda arr: ~select_atom_names(arr, atomname))

    def not_canonical_amino_acids(self):
        return self._update_mask(lambda arr: ~select_canonical_amino_acids(arr))

    def not_canonical_nucleotides(self):
        return self._update_mask(lambda arr: ~select_canonical_nucleotides(arr))

    def not_carbohydrates(self):
        return self._update_mask(lambda arr: ~select_carbohydrates(arr))

    def not_chain_id(self, chain_id):
        return self._update_mask(lambda arr: ~select_chain_id(arr, chain_id))

    def not_element(self, element):
        return self._update_mask(lambda arr: ~select_element(arr, element))

    def not_hetero(self):
        return self._update_mask(lambda arr: ~select_hetero(arr))

    # def not_inscode(self, inscode):
    #     return self._update_mask(lambda arr: ~select_inscode(arr, inscode))

    def not_monoatomic_ions(self):
        return self._update_mask(lambda arr: ~select_monoatomic_ions(arr))

    def not_nucleotides(self):
        return self._update_mask(lambda arr: ~select_nucleotides(arr))

    def not_peptide_backbone(self):
        return self._update_mask(lambda arr: ~select_peptide_backbone(arr))

    def not_phosphate_backbone(self):
        return self._update_mask(lambda arr: ~select_phosphate_backbone(arr))

    def not_polymer(self):
        return self._update_mask(lambda arr: ~select_polymer(arr))

    def not_res_id(self, num):
        return self._update_mask(lambda arr: ~select_res_id(arr, num))

    def not_res_ids(self, nums):
        return self._update_mask(lambda arr: ~select_res_ids(arr, nums))

    def not_res_name(self, res_name):
        return self._update_mask(lambda arr: ~select_res_name(arr, res_name))

    def not_solvent(self):
        return self._update_mask(lambda arr: ~select_solvent(arr))


def select_amino_acids(arr):
    return filter.filter_amino_acids(arr)


def select_atom_names(arr, atom_name):
    return np.isin(arr.get_annotation("atom_name"), atom_name)


def select_canonical_amino_acids(arr):
    return filter.filter_canonical_amino_acids(arr)


def select_canonical_nucleotides(arr):
    return filter.filter_canonical_nucleotides(arr)


def select_carbohydrates(arr):
    return filter.filter_carbohydrates(arr)


def select_chain_id(arr, chain):
    return chain == arr.get_annotation("chain_id")


def select_element(arr, element):
    return element == arr.get_annotation("element")


# def select_first_altloc(arr):
#     return filter.filter_first_altloc(arr)


def select_hetero(arr):
    return arr.get_annotation("hetero")


# def select_highest_occupancy_altloc(arr):
#     return filter.filter_highest_occupancy_altloc(arr)


# def select_inscode(arr, inscode):
#     return inscode == arr.get_annotation("ins_code")


# def select_intersection(arr):
#     return filter.filter_intersection(arr)


def select_ligand(arr):
    """Selects ligand molecules (hetero compounds that are not solvent or ions)"""
    hetero_mask = select_hetero(arr)
    not_solvent_mask = ~select_solvent(arr)
    not_ion_mask = ~select_monoatomic_ions(arr)
    return hetero_mask & not_solvent_mask & not_ion_mask


def select_linear_bond_continuity(arr):
    return filter.filter_linear_bond_continuity(arr)


def select_monoatomic_ions(arr):
    return filter.filter_monoatomic_ions(arr)


def select_nucleotides(arr):
    return filter.filter_nucleotides(arr)


def select_peptide_backbone(arr):
    return filter.filter_peptide_backbone(arr)


def select_phosphate_backbone(arr):
    return filter.filter_phosphate_backbone(arr)


def select_polymer(arr):
    return filter.filter_polymer(arr)


def select_res_id(arr, num):
    return num == arr.get_annotation("res_id")


def select_res_ids(arr, nums):
    return np.isin(arr.get_annotation("res_id"), nums)


def select_res_name(arr, res_name):
    return res_name == arr.get_annotation("res_name")


def select_solvent(arr):
    return filter.filter_solvent(arr)
