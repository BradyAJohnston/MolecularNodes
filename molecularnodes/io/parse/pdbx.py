import biotite.structure.io.pdbx as pdbx
import biotite.structure as struc
import numpy as np

from .molecule import Molecule
from ...assembly.cif import CIFAssemblyParser

class PDBX(Molecule):
    def __init__(self, file_path):
        self.file_path = file_path
        self.file = pdbx.PDBxFile.read(file_path)
        self.structure = self.get_structure()
        self.assemblies = CIFAssemblyParser(self.file).get_assemblies()
        self.n_models = self.structure.shape[0]
        self.n_atoms = self.structure.shape[1]
    
    def get_structure(self):
        array = pdbx.get_structure(self.file, extra_fields = ['b_factor', 'charge', 'occupancy', 'atom_id'])
        array.set_annotation('sec_struct', get_ss_mmcif(array, self.file))
        if not array.bonds:
            array[0].bonds = struc.bonds.connect_via_residue_names(array[0], inter_residue = True)
        
        return array

def ss_id_to_numeric(id: str) -> int:
    "Convert the given ids in the mmmCIF file to 1 AH / 2 BS / 3 Loop integers"
    if "HELX" in id:
        return int(1)
    elif "STRN" in id:
        return int(2)
    else:
        return int(3)

class NoSecondaryStructureError(Exception):
    """Raised when no secondary structure is found"""
    pass

def get_ss_mmcif(mol, file):
    import biotite.structure as struc
    
    conf = file.get_category('struct_conf')
    if not conf:
        raise NoSecondaryStructureError
    starts = conf['beg_auth_seq_id'].astype(int)
    ends = conf['end_auth_seq_id'].astype(int)
    chains = conf['end_auth_asym_id'].astype(str)
    id_label = conf['id'].astype(str)
    
    sheet = file.get_category('struct_sheet_range')
    if sheet:
        starts = np.append(starts, sheet['beg_auth_seq_id'].astype(int))
        ends = np.append(ends, sheet['end_auth_seq_id'].astype(int))
        chains = np.append(chains, sheet['end_auth_asym_id'].astype(str))
        id_label = np.append(id_label, np.repeat('STRN', len(sheet['id'])))
    
    id_int = np.array([ss_id_to_numeric(x) for x in id_label])
    lookup = dict()
    for chain in np.unique(chains):
        arrays = []
        mask = (chain == chains)
        start_sub = starts[mask]
        end_sub = ends[mask]
        id_sub = id_int[mask]
        
        for (start, end, id) in zip(start_sub, end_sub, id_sub):
            idx = np.arange(start, end + 1, dtype = int)
            arr = np.zeros((len(idx), 2), dtype = int)
            arr[:, 0] = idx
            arr[:, 1] = 3
            arr[:, 1] = id
            arrays.append(arr)
        
        lookup[chain] =  dict(np.vstack(arrays).tolist())
    
    ss = []
    
    for i, (chain_id, res_id) in enumerate(zip(mol.chain_id, mol.res_id)):
        ss.append(lookup[chain_id].get(res_id, 3))
    
    arr = np.array(ss, dtype = int)
    arr[~struc.filter_amino_acids(mol)] = 0
    return arr