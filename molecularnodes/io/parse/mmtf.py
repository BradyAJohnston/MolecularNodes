import biotite.structure.io.mmtf as mmtf
import numpy as np

from .assembly import AssemblyParser
from .molecule import Molecule

class MMTF(Molecule):
    def __init__(self, file_path):
        self.file_path = file_path
        self.file = mmtf.MMTFFile.read(self.file_path)
        self.structure = self.get_structure()
        self.assemblies = MMTFAssemblyParser(self.file).get_assemblies()
        self.n_models = self.structure.shape[0]
        self.n_atoms = self.structure.shape[1]

    def get_structure(self):
        array = mmtf.get_structure(
            file = self.file, 
            include_bonds = True, 
            extra_fields = ['b_factor', 'charge', 'occupancy', 'atom_id']
            )
        entity_lookup = get_chain_entity_id(self.file)
        entity_ids = np.array([entity_lookup[x] for x in array.chain_id], dtype = int)
        array.set_annotation('entity_id', entity_ids)
        array.set_annotation('sec_struct', get_secondary_structure(array, self.file))
        return array


def get_secondary_structure(array, file) -> np.array:
    """
    Gets the secondary structure annotation that is included in mmtf files and returns it as a numerical numpy array.

    Parameters:
    -----------
    array : numpy.array
        The molecular coordinates array, from mmtf.get_structure()
    file : mmtf.MMTFFile
        The MMTF file containing the secondary structure information, from mmtf.MMTFFile.read()

    Returns:
    --------
    atom_sse : numpy.array
        Numerical numpy array representing the secondary structure of the molecule.
    
    Description:
    ------------
    This function uses the biotite.structure package to extract the secondary structure information from the MMTF file.
    The resulting secondary structures are `1: Alpha Helix, 2: Beta-sheet, 3: loop`.
    """
    
    from biotite.structure import spread_residue_wise
    
    sec_struct_codes = {
        -1: "X", # undefined
        0 : "I", # pi helix
        1 : "S", # bend
        2 : "H", # alpha helix
        3 : "E", # extended
        4 : "G", # 3-10 helix
        5 : "B", # bridge
        6 : "T", # turn
        7 : "C"  # coil
    }
    
    # convert to 1 AH / 2 BS / 3 LOOP
    dssp_codes_to_int = {
        -1: 0, # undefined
        
        0 : 1, # pi helix
        2 : 1, # alpha helix
        4 : 1, # 3-10 helix
        
        3 : 2, # extended
        5 : 2, # bridge
        
        6 : 3, # turn
        1 : 3, # bend
        7 : 3  # coil
    }
    
    dssp_to_abc = {
        "X" : 0,
        "I" : 3, #"a",
        "G" : 1, #"a",
        "H" : 1, #"a",
        
        "E" : 2, #"b",
        "B" : 2, #"b",
        
        "T" : 3, #"c",
        "S" : 3, #"c",
        "C" : 3  #"c"
    }
    
    try:
        sse = file["secStructList"]
    except KeyError:
        ss_int = np.full(len(array), 3)
        print('Warning: "secStructList" field missing from MMTF file. Defaulting \
            to "loop" for all residues.')
    else:
        pass
        ss_int = np.array(
            [dssp_to_abc.get(sec_struct_codes.get(ss)) for ss in sse], 
            dtype = int
        )
    atom_sse = spread_residue_wise(array, ss_int)
    # atom_sse = spread_residue_wise(array, sse)
    
    return atom_sse

def set_atom_entity_id(mol, file):
    mol.add_annotation('entity_id', int)
    ent_dic = get_chain_entity_id(file)
    
    entity_ids = np.array([ent_dic[x] for x in mol.chain_id])
    
    # entity_ids = chain_entity_id[chain_ids]
    mol.set_annotation('entity_id', entity_ids)
    return entity_ids

def get_chain_entity_id(file):
    entities = file['entityList']
    chain_names = file['chainNameList']    
    ent_dic = {}
    for i, ent in enumerate(entities):
        for chain_idx in ent['chainIndexList']:
            chain_id = chain_names[chain_idx]
            if  chain_id in ent_dic.keys():
                next
            else:
                ent_dic[chain_id] = i
    
    return ent_dic

class MMTFAssemblyParser(AssemblyParser):
    ### Implementation adapted from ``biotite.structure.io.mmtf.assembly``

    def __init__(self, mmtf_file):
        self._file = mmtf_file
    

    def list_assemblies(self):
        import biotite.structure.io.mmtf as mmtf
        return mmtf.list_assemblies(self._file)
    

    def get_transformations(self, assembly_id):
        import biotite
        # Find desired assembly
        selected_assembly = None
        if "bioAssemblyList" not in self._file:
            raise biotite.InvalidFileError(
                "File does not contain assembly information "
                "(missing 'bioAssemblyList')"
            )
        for assembly in self._file["bioAssemblyList"]:
            current_assembly_id = assembly["name"]
            transform_list = assembly["transformList"]
            if current_assembly_id == assembly_id:
                selected_assembly = transform_list
                break
        if selected_assembly is None:
            raise KeyError(
                f"The assembly ID '{assembly_id}' is not found"
            )

        # Parse transformations from assembly
        transformations = []
        for transform in selected_assembly:
            matrix = np.array(transform["matrix"]).reshape(4, 4)
            chain_ids = np.array(self._file["chainNameList"], dtype="U4")
            affected_chain_ids = chain_ids[transform["chainIndexList"]]
            
            transformations.append((
                affected_chain_ids.tolist(), 
                matrix.tolist()
            ))

        return transformations
    
    def get_assemblies(self):
        assembly_dict = {}
        for assembly_id in self.list_assemblies():
            assembly_dict[assembly_id] = self.get_transformations(assembly_id)
        
        return assembly_dict