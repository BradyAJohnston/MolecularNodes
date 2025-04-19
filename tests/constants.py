import os
from pathlib import Path

codes = ["4ozs", "8H1B", "1BNA", "8U8W"]
attributes = [
    "b_factor",
    "occupancy",
    "vdw_radii",
    "lipophobicity",
    "charge",
    "res_id",
    "res_name",
    "atomic_number",
    "chain_id",
    "entity_id",
    "atom_id",
    "atom_name",
    "sec_struct",
    "Color",
    "position",
    "is_backbone",
    "is_side_chain",
    "is_alpha_carbon",
    "is_solvent",
    "is_nucleic",
    "is_peptide",
    "is_hetero",
    "is_carb",
    "bond_type",
    "mass",
    "ures_id",
]
data_dir = Path(os.path.abspath(Path(__file__).parent / "data"))
