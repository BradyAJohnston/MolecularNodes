import numpy as np
from biotite.structure import filter, AtomArray, AtomArrayStack
from typing import Callable


def select_amino_acids(arr):
    return filter.filter_amino_acids(arr)


def select_atom_names(arr, atom_name):
    if isinstance(atom_name, str):
        atom_name = [atom_name]
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


def select_res_id(arr, nums):
    return np.isin(arr.get_annotation("res_id"), nums)


def select_res_name(arr, res_name):
    if isinstance(res_name, str):
        res_name = [res_name]
    return np.isin(arr.get_annotation("res_name"), res_name)


def select_solvent(arr):
    return filter.filter_solvent(arr)
