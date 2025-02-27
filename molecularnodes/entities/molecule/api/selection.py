import numpy as np
from biotite.structure import filter
from types import MethodType


class StructureSelector:
    def _update_mask(self, new_mask):
        self.mask = new_mask if self.mask is None else (self.mask & new_mask)
        return self

    def amino_acids(self):
        return self._update_mask(select_amino_acids(self.array))

    def atomname(self, atomname):
        return self._update_mask(select_atomname(self.array, atomname))

    def canonical_amino_acids(self):
        return self._update_mask(select_canonical_amino_acids(self.array))

    def canonical_nucleotides(self):
        return self._update_mask(select_canonical_nucleotides(self.array))

    def carbohydrates(self):
        return self._update_mask(select_carbohydrates(self.array))

    def chain(self, chain_id):
        return self._update_mask(select_chain(self.array, chain_id))

    def element(self, element):
        return self._update_mask(select_element(self.array, element))

    def first_altloc(self):
        return self._update_mask(select_first_altloc(self.array))

    def get_mask(self):
        return self.mask

    def get_selection(self):
        """Returns the structure filtered by the current mask"""
        if self.mask is None:
            return self.array
        return self.array[self.mask]

    def hetero(self):
        return self._update_mask(select_hetero(self.array))

    def highest_occupancy_altloc(self):
        return self._update_mask(select_highest_occupancy_altloc(self.array))

    def inscode(self, inscode):
        return self._update_mask(select_inscode(self.array, inscode))

    def intersection(self):
        return self._update_mask(select_intersection(self.array))

    def ligand(self):
        return self._update_mask(select_ligand(self.array))

    def linear_bond_continuity(self):
        return self._update_mask(select_linear_bond_continuity(self.array))

    def monoatomic_ions(self):
        return self._update_mask(select_monoatomic_ions(self.array))

    def nucleotides(self):
        return self._update_mask(select_nucleotides(self.array))

    def peptide_backbone(self):
        return self._update_mask(select_peptide_backbone(self.array))

    def phosphate_backbone(self):
        return self._update_mask(select_phosphate_backbone(self.array))

    def polymer(self):
        return self._update_mask(select_polymer(self.array))

    def resid(self, num):
        return self._update_mask(select_resid(self.array, num))

    def resids(self, nums):
        return self._update_mask(select_resids(self.array, nums))

    def resname(self, res_name):
        return self._update_mask(select_resname(self.array, res_name))

    def solvent(self):
        return self._update_mask(select_solvent(self.array))

    def not_amino_acids(self):
            return self._update_mask(~select_amino_acids(self.array))

    def not_atomname(self, atomname):
        return self._update_mask(~select_atomname(self.array, atomname))

    def not_canonical_amino_acids(self):
        return self._update_mask(~select_canonical_amino_acids(self.array))

    def not_canonical_nucleotides(self):
        return self._update_mask(~select_canonical_nucleotides(self.array))

    def not_carbohydrates(self):
        return self._update_mask(~select_carbohydrates(self.array))

    def not_chain(self, chain_id):
        return self._update_mask(~select_chain(self.array, chain_id))

    def not_element(self, element):
        return self._update_mask(~select_element(self.array, element))

    def not_hetero(self):
        return self._update_mask(~select_hetero(self.array))

    def not_inscode(self, inscode):
        return self._update_mask(~select_inscode(self.array, inscode))

    def not_monoatomic_ions(self):
        return self._update_mask(~select_monoatomic_ions(self.array))

    def not_nucleotides(self):
        return self._update_mask(~select_nucleotides(self.array))

    def not_peptide_backbone(self):
        return self._update_mask(~select_peptide_backbone(self.array))

    def not_phosphate_backbone(self):
        return self._update_mask(~select_phosphate_backbone(self.array))

    def not_polymer(self):
        return self._update_mask(~select_polymer(self.array))

    def not_resid(self, num):
        return self._update_mask(~select_resid(self.array, num))

    def not_resids(self, nums):
        return self._update_mask(~select_resids(self.array, nums))

    def not_resname(self, res_name):
        return self._update_mask(~select_resname(self.array, res_name))

    def not_solvent(self):
        return self._update_mask(~select_solvent(self.array))


def select_amino_acids(arr):
    return filter.filter_amino_acids(arr)

def select_atomname(arr, atomname):
    return atomname == arr.get_annotation("atom_name")

def select_canonical_amino_acids(arr):
    return filter.filter_canonical_amino_acids(arr)

def select_canonical_nucleotides(arr):
    return filter.filter_canonical_nucleotides(arr)

def select_carbohydrates(arr):
    return filter.filter_carbohydrates(arr)

def select_chain(arr, chain):
    return chain == arr.get_annotation("chain_id")

def select_element(arr, element):
    return element == arr.get_annotation("element")

def select_first_altloc(arr):
    return filter.filter_first_altloc(arr)

def select_hetero(arr):
    return True == arr.get_annotation("hetero")

def select_highest_occupancy_altloc(arr):
    return filter.filter_highest_occupancy_altloc(arr)

def select_inscode(arr, inscode):
    return inscode == arr.get_annotation("ins_code")

def select_intersection(arr):
    return filter.filter_intersection(arr)

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

def select_resid(arr, num):
    return num == arr.get_annotation("res_id")

def select_resids(arr, nums):
    return np.isin(arr.get_annotation("res_id"), nums)

def select_resname(arr, res_name):
    return res_name == arr.get_annotation("res_name")

def select_solvent(arr):
    return filter.filter_solvent(arr)
