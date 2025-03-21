import numpy as np
from biotite.structure import filter, AtomArray, AtomArrayStack
from typing import Callable


class Selector:
    """
    A helper to create selections for Molecules and AtomArrays.

    The selection (self.mask) is not computed or returned until the `evaluate_on_array`
    method is called. Until then methods are stored for later evaluation.

    Parameters
    ----------
    None

    Attributes
    ----------
    mask : ndarray or None
        Boolean array for the selection on the most recently evaluated array.
    pending_selections : list
        List of selection operations to be applied once `evaluate_on_array` is called.
    """

    def __init__(self):
        self.mask: np.ndarray | None = None
        self.pending_selections: list[Callable] = []

    def _update_mask(self, operation):
        """
        Add a selection operation to the pending list.

        Parameters
        ----------
        operation : callable
            Selection operation to add

        Returns
        -------
        self : Selector
            Returns self for method chaining
        """
        self.pending_selections.append(operation)
        return self

    def reset(self):
        """
        Reset all pending selections and the mask

        Returns
        -------
        self : Selector
            Returns self for method chaining
        """
        self.pending_selections = []
        self.mask = None
        return self

    def evaluate_on_array(self, array: AtomArray | AtomArrayStack) -> np.ndarray:
        """Evaluate this selection on the AtomArray.

        Parameters
        ----------
        array : AtomArray or AtomArrayStack
            The atomic structure to evaluate the selection on.

        Returns
        -------
        ndarray
            Boolean mask array indicating which atoms match the selection criteria.

        Notes
        -----
        All of the selection operations that have been staged for this Selector are
        evaluated and combined with a logical AND, using the AtomArray as input.
        """
        if isinstance(array, AtomArrayStack):
            array = array[0]  # type: ignore

        self.mask = np.ones(array.array_length(), dtype=bool)
        if not self.pending_selections:
            return self.mask

        for operation in self.pending_selections:
            self.mask &= operation(array)
        return self.mask  # type: ignore

    def amino_acids(self):
        """Select amino acid residues.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(select_amino_acids)

    def atom_name(self, atom_name: str | list[str] | tuple[str, ...] | np.ndarray):
        """Select atoms by their name.

        Parameters
        ----------
        atom_name : str or list of str or tuple of str or ndarray
            The atom name(s) to select.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: select_atom_names(arr, atom_name))

    def is_canonical_amino_acid(self):
        """Select canonical amino acid residues.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(select_canonical_amino_acids)

    def is_canonical_nucleotide(self):
        """Select canonical nucleotide residues.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(select_canonical_nucleotides)

    def is_carbohydrate(self):
        """Select carbohydrate residues.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(select_carbohydrates)

    def chain_id(self, chain_id: list[str] | tuple[str, ...] | np.ndarray):
        """Select atoms by chain identifier.

        Parameters
        ----------
        chain_id : list of str or tuple of str or ndarray
            The chain identifier(s) to select.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: select_chain_id(arr, chain_id))

    def element(self, element: list[str] | tuple[str, ...] | np.ndarray):
        """Select atoms by element symbol.

        Parameters
        ----------
        element : list of str or tuple of str or ndarray
            The element symbol(s) to select.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: select_element(arr, element))

    def is_hetero(self):
        """Select hetero atoms.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(select_hetero)

    def is_ligand(self):
        """Select ligand atoms.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(select_ligand)

    def linear_bond_continuity(self):
        """Select atoms with linear bond continuity.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(select_linear_bond_continuity)

    def is_monoatomic_ion(self):
        """Select monoatomic ions.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(select_monoatomic_ions)

    def is_nucleotide(self):
        """Select nucleotide residues.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(select_nucleotides)

    def peptide_backbone(self):
        """Select peptide backbone atoms.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(select_peptide_backbone)

    def phosphate_backbone(self):
        """Select phosphate backbone atoms.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(select_phosphate_backbone)

    def polymer(self):
        """Select polymer atoms.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(select_polymer)

    def res_id(self, num):
        """Select atoms by residue ID.

        Parameters
        ----------
        num : int or list of int
            The residue ID(s) to select.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: select_res_id(arr, num))

    def res_name(self, res_name):
        """Select atoms by residue name.

        Parameters
        ----------
        res_name : str or list of str
            The residue name(s) to select.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: select_res_name(arr, res_name))

    def solvent(self):
        """Select solvent atoms.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(select_solvent)

    def not_amino_acids(self):
        """Select non-amino acid residues.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: ~select_amino_acids(arr))

    def not_atom_names(self, atomname):
        """Select atoms not matching the specified atom names.

        Parameters
        ----------
        atomname : str or list of str
            The atom name(s) to exclude.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: ~select_atom_names(arr, atomname))

    def not_canonical_amino_acids(self):
        """Select non-canonical amino acid residues.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: ~select_canonical_amino_acids(arr))

    def not_canonical_nucleotides(self):
        """Select non-canonical nucleotide residues.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: ~select_canonical_nucleotides(arr))

    def not_carbohydrates(self):
        """Select non-carbohydrate residues.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: ~select_carbohydrates(arr))

    def not_chain_id(self, chain_id):
        """Select atoms not in the specified chains.

        Parameters
        ----------
        chain_id : str or list of str
            The chain identifier(s) to exclude.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: ~select_chain_id(arr, chain_id))

    def not_element(self, element):
        """Select atoms not matching the specified elements.

        Parameters
        ----------
        element : str or list of str
            The element symbol(s) to exclude.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: ~select_element(arr, element))

    def not_hetero(self):
        """Select non-hetero atoms.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: ~select_hetero(arr))

    def not_monoatomic_ions(self):
        """Select non-monoatomic ion atoms.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: ~select_monoatomic_ions(arr))

    def not_nucleotides(self):
        """Select non-nucleotide residues.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: ~select_nucleotides(arr))

    def not_peptide_backbone(self):
        """Select non-peptide backbone atoms.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: ~select_peptide_backbone(arr))

    def not_phosphate_backbone(self):
        """Select non-phosphate backbone atoms.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: ~select_phosphate_backbone(arr))

    def not_polymer(self):
        """Select non-polymer atoms.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: ~select_polymer(arr))

    def not_res_id(self, num):
        """Select atoms not matching the specified residue IDs.

        Parameters
        ----------
        num : int or list of int
            The residue ID(s) to exclude.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: ~select_res_id(arr, num))

    def not_res_name(self, res_name):
        """Select atoms not matching the specified residue names.

        Parameters
        ----------
        res_name : str or list of str
            The residue name(s) to exclude.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: ~select_res_name(arr, res_name))

    def not_solvent(self):
        """Select non-solvent atoms.

        Returns
        -------
        Selector
            Returns self for method chaining.
        """
        return self._update_mask(lambda arr: ~select_solvent(arr))


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
