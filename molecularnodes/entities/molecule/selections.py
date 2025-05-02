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


def select_peptide(arr: AtomArray | AtomArrayStack) -> np.ndarray:
    return np.logical_or(
        select_amino_acids(arr),
        select_canonical_amino_acids(arr),
    )


def select_carbohydrates(arr: AtomArray | AtomArrayStack) -> np.ndarray:
    return filter.filter_carbohydrates(arr)


def select_chain_id(
    arr: AtomArray | AtomArrayStack,
    chain: str | list[str] | tuple[str, ...] | np.ndarray,
) -> np.ndarray:
    return np.isin(arr.get_annotation("chain_id"), chain)


def select_element(
    arr: AtomArray | AtomArrayStack,
    element: str | list[str] | tuple[str, ...] | np.ndarray,
) -> np.ndarray:
    return np.isin(arr.get_annotation("element"), element)


def select_hetero(arr: AtomArray | AtomArrayStack) -> np.ndarray:
    return arr.get_annotation("hetero")


def select_ligand(arr: AtomArray | AtomArrayStack) -> np.ndarray:
    """Selects ligand molecules (hetero compounds that are not solvent or ions)"""
    hetero_mask = select_hetero(arr)
    not_solvent_mask = ~select_solvent(arr)
    not_ion_mask = ~select_monoatomic_ions(arr)
    return hetero_mask & not_solvent_mask & not_ion_mask


def select_linear_bond_continuity(arr: AtomArray | AtomArrayStack) -> np.ndarray:
    return filter.filter_linear_bond_continuity(arr)


def select_monoatomic_ions(arr: AtomArray | AtomArrayStack) -> np.ndarray:
    return filter.filter_monoatomic_ions(arr)


def select_nucleotides(arr: AtomArray | AtomArrayStack) -> np.ndarray:
    return filter.filter_nucleotides(arr)


def select_peptide_backbone(arr: AtomArray | AtomArrayStack) -> np.ndarray:
    return filter.filter_peptide_backbone(arr)


def select_phosphate_backbone(arr: AtomArray | AtomArrayStack) -> np.ndarray:
    return filter.filter_phosphate_backbone(arr)


def select_alpha_carbon(arr: AtomArray | AtomArrayStack) -> np.ndarray:
    return arr.get_annotation("atom_name") == "CA"


def select_backbone(arr: AtomArray | AtomArrayStack) -> np.ndarray:
    """
    Get the atoms that appear in peptide backbone or nucleic acid phosphate backbones.
    Filter differs from the Biotite's `struc.filter_peptide_backbone()` in that this
    includes the peptide backbone oxygen atom, which biotite excludes. Additionally
    this selection also includes all of the atoms from the ribose in nucleic acids,
    and the other phosphate oxygens.
    """
    backbone_atom_names = [
        # Peptide backbone atoms
        "N",
        "C",
        "CA",
        "H",
        "HA",
        "O",
        # Continuous nucleic backbone atoms
        "P",
        "O5'",
        "C5'",
        "C4'",
        "C3'",
        "O3'",
        # Alternative names for phosphate O's
        "O1P",
        "OP1",
        "O2P",
        "OP2",
        # Remaining ribose atoms
        "O4'",
        "C1'",
        "C2'",
        "O2'",
    ]

    is_backbone_atom = np.isin(arr.get_annotation("atom_name"), backbone_atom_names)
    is_not_solvent = np.logical_not(filter.filter_solvent(arr))

    return np.logical_and(is_backbone_atom, is_not_solvent)


def select_polymer(arr: AtomArray | AtomArrayStack) -> np.ndarray:
    return filter.filter_polymer(arr)


def select_side_chain(arr: AtomArray | AtomArrayStack) -> np.ndarray:
    not_backbone = ~select_backbone(arr)
    is_polymer = np.logical_or(select_nucleotides(arr), select_peptide(arr))
    # is_ca = select_alpha_carbon(arr)
    final = np.logical_and(not_backbone, is_polymer)
    # TODO: to match with previous selection atributes we don't include the CA atom, but
    # this will be changed in a future PR
    # final = np.logical_or(final, is_ca)
    return final


def select_res_id(
    arr: AtomArray | AtomArrayStack,
    nums: int | range | list[int] | tuple[int, ...] | np.ndarray,
):
    return np.isin(arr.get_annotation("res_id"), nums)


def select_res_name(
    arr: AtomArray | AtomArrayStack,
    res_name: str | list[str] | tuple[str, ...] | np.ndarray,
) -> np.ndarray:
    return np.isin(arr.get_annotation("res_name"), res_name)


def select_solvent(arr: AtomArray | AtomArrayStack) -> np.ndarray:
    return filter.filter_solvent(arr)
