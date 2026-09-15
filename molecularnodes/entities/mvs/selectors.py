"""
Translate MolViewSpec component selectors into boolean atom masks.

MVS selects atoms either with a static selector (``"polymer"``, ``"ligand"``,
...) or with one or more ``ComponentExpression`` objects whose fields are
AND-ed together (a list of expressions is OR-ed). Masks are computed here in
Python against the entity's universe and import-time attributes, and the
caller stores them as named boolean attributes for the style branches to use.
"""

from typing import TYPE_CHECKING, Any
import numpy as np

if TYPE_CHECKING:
    from ..molecule.base import Molecule


def _attribute(mol: "Molecule", name: str) -> np.ndarray | None:
    try:
        return mol.named_attribute(name)
    except Exception:
        return None


def _residue_sizes(mol: "Molecule") -> np.ndarray:
    """Per-atom size of the residue each atom belongs to."""
    atoms = mol.universe.atoms
    sizes = np.array([residue.atoms.n_atoms for residue in atoms.residues])
    return sizes[atoms.resindices - atoms.resindices.min()]


def static_selector_mask(
    mol: "Molecule", selector: str, warnings_out: list[str]
) -> np.ndarray:
    """Mask for one of the enumerated MVS selectors."""
    n_atoms = mol.universe.atoms.n_atoms

    def attr(name: str) -> np.ndarray:
        values = _attribute(mol, name)
        if values is None:
            return np.zeros(n_atoms, dtype=bool)
        return values.astype(bool)

    match selector:
        case "all":
            return np.ones(n_atoms, dtype=bool)
        case "polymer":
            return attr("is_peptide") | attr("is_nucleic")
        case "protein":
            return attr("is_peptide")
        case "nucleic":
            return attr("is_nucleic")
        case "branched":
            return attr("is_carb")
        case "water":
            return attr("is_solvent")
        case "ion":
            # single-atom hetero residues that are not water
            return attr("is_hetero") & ~attr("is_solvent") & (_residue_sizes(mol) == 1)
        case "ligand":
            return (
                attr("is_hetero")
                & ~attr("is_solvent")
                & ~attr("is_carb")
                & (_residue_sizes(mol) != 1)
            )
        case "coarse":
            warnings_out.append(
                "the 'coarse' selector is not supported; selecting nothing"
            )
            return np.zeros(n_atoms, dtype=bool)
        case _:
            warnings_out.append(f"unknown selector '{selector}'; selecting nothing")
            return np.zeros(n_atoms, dtype=bool)


def expression_mask(
    mol: "Molecule",
    expression: dict[str, Any],
    label_fields: dict[str, np.ndarray],
    warnings_out: list[str],
) -> np.ndarray:
    """
    Mask for a single ``ComponentExpression``: all given fields are AND-ed.

    ``label_*`` identifier fields are resolved from ``label_fields`` (a
    secondary parse of the source file) when available, and otherwise fall
    back to the corresponding ``auth_*`` values with a warning.
    """
    atoms = mol.universe.atoms
    mask = np.ones(atoms.n_atoms, dtype=bool)

    def label_or_auth(field: str, auth_values: np.ndarray) -> np.ndarray:
        if field in label_fields:
            return label_fields[field]
        warnings_out.append(
            f"'{field}' is not recorded for this structure; using the "
            "author-assigned equivalent instead"
        )
        return auth_values

    for field, value in expression.items():
        if value is None:
            continue
        match field:
            case "auth_asym_id":
                mask &= atoms.chainIDs == value
            case "label_asym_id":
                mask &= label_or_auth("label_asym_id", atoms.chainIDs) == str(value)
            case "label_entity_id":
                entity_ids = list(mol.props.entity_ids)
                entity_attr = _attribute(mol, "entity_id")
                if entity_attr is None or str(value) not in entity_ids:
                    mask &= False
                else:
                    mask &= entity_attr == entity_ids.index(str(value))
            case "auth_seq_id":
                mask &= atoms.resids == value
            case "label_seq_id":
                mask &= label_or_auth("label_seq_id", atoms.resids) == value
            case "beg_auth_seq_id":
                mask &= atoms.resids >= value
            case "end_auth_seq_id":
                mask &= atoms.resids <= value
            case "beg_label_seq_id":
                mask &= label_or_auth("label_seq_id", atoms.resids) >= value
            case "end_label_seq_id":
                mask &= label_or_auth("label_seq_id", atoms.resids) <= value
            case "label_comp_id" | "auth_comp_id":
                mask &= atoms.resnames == value
            case "label_atom_id" | "auth_atom_id":
                mask &= atoms.names == value
            case "type_symbol":
                mask &= np.char.upper(atoms.elements.astype(str)) == str(value).upper()
            case "atom_id":
                atom_ids = _attribute(mol, "atom_id")
                if atom_ids is None:
                    mask &= False
                else:
                    mask &= atom_ids == value
            case "atom_index":
                index_mask = np.zeros(atoms.n_atoms, dtype=bool)
                if 0 <= value < atoms.n_atoms:
                    index_mask[value] = True
                mask &= index_mask
            case "residue_index":
                ures_id = _attribute(mol, "ures_id")
                if ures_id is None:
                    mask &= False
                else:
                    mask &= ures_id == value
            case "pdbx_PDB_ins_code":
                icodes = getattr(atoms, "icodes", None)
                if icodes is None:
                    warnings_out.append(
                        "insertion codes are not recorded for this structure; "
                        "ignoring 'pdbx_PDB_ins_code'"
                    )
                else:
                    mask &= icodes == value
            case _:
                warnings_out.append(
                    f"unsupported component expression field '{field}'; ignoring it"
                )

    return mask


def component_mask(
    mol: "Molecule",
    selector: Any,
    label_fields: dict[str, np.ndarray],
    warnings_out: list[str],
) -> np.ndarray:
    """
    Mask for a component's ``selector`` parameter.

    A string is a static selector, a dict is a single expression, and a list
    of dicts is OR-ed together. ``None`` selects everything (the MVS default).
    """
    if selector is None:
        selector = "all"
    if isinstance(selector, str):
        return static_selector_mask(mol, selector, warnings_out)
    if isinstance(selector, dict):
        return expression_mask(mol, selector, label_fields, warnings_out)
    if isinstance(selector, (list, tuple)):
        mask = np.zeros(mol.universe.atoms.n_atoms, dtype=bool)
        for expression in selector:
            mask |= component_mask(mol, expression, label_fields, warnings_out)
        return mask
    warnings_out.append(f"unsupported selector {selector!r}; selecting nothing")
    return np.zeros(mol.universe.atoms.n_atoms, dtype=bool)
