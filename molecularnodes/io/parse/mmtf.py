import numpy as np
import warnings

from .assembly import AssemblyParser
from .molecule import Molecule


class MMTF(Molecule):
    def __init__(self, file_path):
        super().__init__()
        self.file_path = file_path
        self.file = self._read()
        self.array = self._get_structure()
        self.n_models = self.array.shape[0]
        self.n_atoms = self.array.shape[1]
        self.entity_ids = self._entity_ids()
        self.chain_ids = self._chain_ids()

    def _read(self):
        import biotite.structure.io.mmtf as mmtf
        return mmtf.MMTFFile.read(self.file_path)

    def _get_structure(self):
        import biotite.structure.io.mmtf as mmtf
        array = mmtf.get_structure(
            file=self.file,
            include_bonds=True,
            extra_fields=['b_factor', 'charge', 'occupancy', 'atom_id']
        )

        array.set_annotation('entity_id', _get_entity_id(array, self.file))
        array.set_annotation(
            'sec_struct', _get_secondary_structure(array, self.file))
        return array

    def _assemblies(self):
        return MMTFAssemblyParser(self.file).get_assemblies()

    def _entity_ids(self):
        try:
            return [item['description'] for item in self.file['entityList']]
        except KeyError:
            return None


def _get_secondary_structure(array, file) -> np.array:
    """
    Gets the secondary structure annotation that is included in mmtf files and returns it as a numerical numpy array.

    Parameters:
    -----------
    array : numpy.array
        The AtomArray from mmtf.get_structure(mmtf.MMTFFile).
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
        -1: "X",  # undefined
        0: "I",  # pi helix
        1: "S",  # bend
        2: "H",  # alpha helix
        3: "E",  # extended
        4: "G",  # 3-10 helix
        5: "B",  # bridge
        6: "T",  # turn
        7: "C"  # coil
    }

    code_to_abc_int = {
        "X": 0,
        "I": 3,  # "a",
        "G": 1,  # "a",
        "H": 1,  # "a",

        "E": 2,  # "b",
        "B": 2,  # "b",

        "T": 3,  # "c",
        "S": 3,  # "c",
        "C": 3  # "c"
    }

    try:
        secondary_structure = file["secStructList"]
    except KeyError:
        ss_int = np.full(len(array), 3)
        print('Warning: "secStructList" field missing from MMTF file. Defaulting \
            to "loop" for all residues.')
    else:
        ss_str = [sec_struct_codes[x] for x in secondary_structure]
        ss_int = np.array([code_to_abc_int[c] for c in ss_str])

    return spread_residue_wise(array, ss_int)


def _get_entity_id(array, file):
    try:
        entities = file['entityList']
    except KeyError:
        warnings.warn("No entity information in the file.")
        return np.full(len(array), 0)
    names = file['chainNameList']
    entity_lookup = {}
    for i, entity in enumerate(entities):
        for j in entity['chainIndexList']:
            chain_id = names[j]
            if chain_id not in entity_lookup.keys():
                entity_lookup[chain_id] = i

    return np.array([entity_lookup[chain_id] for chain_id in array.chain_id])


class MMTFAssemblyParser(AssemblyParser):
    # Implementation adapted from ``biotite.structure.io.mmtf.assembly``

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
