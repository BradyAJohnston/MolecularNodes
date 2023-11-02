from pathlib import Path
import os

codes = ['4ozs', '8H1B', '1BNA', '8U8W']
attributes = ['b_factor', 'vdw_radii', 'lipophobicity', 'charge', 'res_id', 'res_name', 'atomic_number', 'chain_id', 'entity_id', 'atom_name', 'sec_struct', 'Color', 'position', '.select_vert', 'is_backbone', 'is_alpha_carbon', 'is_solvent', 'is_nucleic', 'is_peptide', 'is_hetero', 'is_carb', 'bond_type', '.edge_verts', '.select_edge', 'sharp_face']
test_data_directory = Path(os.path.abspath(Path(__file__).parent / 'data'))
