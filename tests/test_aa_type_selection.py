"""Test amino acid type selection data layer (non-Blender tests)."""

import sys
from pathlib import Path

# Add assets directory to path to import data.py without triggering bpy import
sys.path.insert(0, str(Path(__file__).parent.parent / "molecularnodes" / "assets"))
import data


def test_residue_type_data_consistency():
    """Verify that all standard amino acids have a res_type defined."""
    standard_aa = [
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    ]

    for aa in standard_aa:
        assert aa in data.residues, f"{aa} missing from residues dict"
        assert "res_type" in data.residues[aa], f"{aa} missing res_type"
        assert data.residues[aa]["res_type"] in [
            "polar", "apolar", "acid", "basic", "aromatic"
        ], f"{aa} has invalid res_type: {data.residues[aa]['res_type']}"


def test_residue_type_classification():
    """Verify residue type assignments match biochemical classifications."""
    polar = ["ASN", "CYS", "GLN", "HIS", "SER", "THR"]
    apolar = ["ALA", "GLY", "ILE", "LEU", "MET", "PRO", "VAL"]
    acidic = ["ASP", "GLU"]
    basic = ["ARG", "LYS"]
    aromatic = ["PHE", "TRP", "TYR"]

    for aa in polar:
        assert data.residues[aa]["res_type"] == "polar", f"{aa} should be polar"

    for aa in apolar:
        assert data.residues[aa]["res_type"] == "apolar", f"{aa} should be apolar"

    for aa in acidic:
        assert data.residues[aa]["res_type"] == "acid", f"{aa} should be acid"

    for aa in basic:
        assert data.residues[aa]["res_type"] == "basic", f"{aa} should be basic"

    for aa in aromatic:
        assert data.residues[aa]["res_type"] == "aromatic", f"{aa} should be aromatic"


def test_residue_name_to_number_mapping():
    """Verify that res_name_num is correctly assigned and unique for standard AAs."""
    standard_aa = [
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    ]

    res_nums = set()
    for aa in standard_aa:
        res_num = data.residues[aa]["res_name_num"]
        assert res_num not in res_nums, f"Duplicate res_name_num {res_num} for {aa}"
        res_nums.add(res_num)
        assert 0 <= res_num <= 19, f"{aa} res_name_num out of range: {res_num}"


def test_select_aa_type_index_mapping():
    """Verify that Index Switch mappings align with res_name_num for each type."""
    # Expected mappings based on res_name_num from data.residues
    polar_indices = {2, 4, 6, 8, 15, 16}  # ASN, CYS, GLN, HIS, SER, THR
    apolar_indices = {0, 7, 9, 10, 12, 14, 19}  # ALA, GLY, ILE, LEU, MET, PRO, VAL
    acidic_indices = {3, 5}  # ASP, GLU
    basic_indices = {1, 11}  # ARG, LYS
    aromatic_indices = {13, 17, 18}  # PHE, TRP, TYR

    # Verify against actual data
    for aa, info in data.residues.items():
        if aa not in ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY",
                      "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
                      "THR", "TRP", "TYR", "VAL"]:
            continue

        res_num = info["res_name_num"]
        res_type = info["res_type"]

        if res_type == "polar":
            assert res_num in polar_indices, f"{aa} (polar) has res_num {res_num} not in polar_indices"
        elif res_type == "apolar":
            assert res_num in apolar_indices, f"{aa} (apolar) has res_num {res_num} not in apolar_indices"
        elif res_type == "acid":
            assert res_num in acidic_indices, f"{aa} (acid) has res_num {res_num} not in acidic_indices"
        elif res_type == "basic":
            assert res_num in basic_indices, f"{aa} (basic) has res_num {res_num} not in basic_indices"
        elif res_type == "aromatic":
            assert res_num in aromatic_indices, f"{aa} (aromatic) has res_num {res_num} not in aromatic_indices"


if __name__ == "__main__":
    test_residue_type_data_consistency()
    print("✓ test_residue_type_data_consistency passed")

    test_residue_type_classification()
    print("✓ test_residue_type_classification passed")

    test_residue_name_to_number_mapping()
    print("✓ test_residue_name_to_number_mapping passed")

    test_select_aa_type_index_mapping()
    print("✓ test_select_aa_type_index_mapping passed")

    print("\nAll tests passed!")
