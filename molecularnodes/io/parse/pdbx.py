import numpy as np
import itertools

from .molecule import Molecule
from ... import utils
from .assembly import AssemblyParser

class PDBX(Molecule):
    def __init__(self, file_path, extra_fields = None, sec_struct=True):
        self.file_path = file_path
        self.file = self._read()
        self.structure = self._get_structure(extra_fields=extra_fields, sec_struct=sec_struct)
        self.n_models = self._n_models()
        self.n_atoms = self._n_atoms()
    def _read(self):
        import biotite.structure.io.pdbx as pdbx
        return pdbx.PDBxFile.read(self.file_path)
    
    def _get_structure(self, extra_fields: str = None, sec_struct=True):
        import biotite.structure.io.pdbx as pdbx
        import biotite.structure as struc
        from biotite import InvalidFileError
        fields = ['b_factor', 'charge', 'occupancy', 'atom_id']
        if extra_fields:
            [fields.append(x) for x in extra_fields]
        
        # if the 'atom_site' doesn't exist then it will just be a small molecule
        # which can be extracted with the get_component()
        try:
            array = pdbx.get_structure(self.file, extra_fields = ['b_factor', 'charge', 'occupancy', 'atom_id'])
            try:
                array.set_annotation('sec_struct', get_ss_mmcif(array, self.file))
            except NoSecondaryStructureError:
                pass
        except InvalidFileError:
            array = pdbx.get_component(self.file)
        
        if not array.bonds:
            array[0].bonds = struc.bonds.connect_via_residue_names(array[0], inter_residue = True)
        
        return array
    
    def _n_models(self):
        import biotite.structure as struc
        if isinstance(self.structure, struc.AtomArray):
            return 1
        else:
            self.structure.shape[0]
    
    def _n_atoms(self):
        import biotite.structure as struc
        arr = self.structure
        if isinstance(self.structure, struc.AtomArray):
            return arr.shape[0]
        else:
            return arr.shape[1]

    def _assemblies(self):
        return CIFAssemblyParser(self.file).get_assemblies()

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
    
    def get_ss(chain_id, res_id):
        try:
            return lookup[chain_id].get(res_id, 3)
        except KeyError:
            return 0
    
    arr = np.array([
            get_ss(chain_id, res_id) for chain_id, res_id in zip(mol.chain_id, mol.res_id)
        ], dtype = int)
    
    arr[~struc.filter_amino_acids(mol)] = 0
    return arr

class CIFAssemblyParser(AssemblyParser):
    ### Implementation adapted from ``biotite.structure.io.pdbx.convert``

    def __init__(self, mmtf_file):
        self._file = mmtf_file
    

    def list_assemblies(self):
        import biotite.structure.io.pdbx as pdbx    
        return list(pdbx.list_assemblies(self._file).keys())
    

    def get_transformations(self, assembly_id):
        import biotite
        assembly_gen_category = self._file.get_category(
            "pdbx_struct_assembly_gen", expect_looped=True
        )
        if assembly_gen_category is None:
            raise biotite.InvalidFileError(
                "File has no 'pdbx_struct_assembly_gen' category"
            )

        struct_oper_category = self._file.get_category(
            "pdbx_struct_oper_list", expect_looped=True
        )
        if struct_oper_category is None:
            raise biotite.InvalidFileError(
                "File has no 'pdbx_struct_oper_list' category"
            )

        if assembly_id not in assembly_gen_category["assembly_id"]:
            raise KeyError(f"File has no Assembly ID '{assembly_id}'")

        # Extract all possible transformations indexed by operation ID
        transformation_dict = _get_transformations(struct_oper_category)

        # Get necessary transformations and the affected chain IDs
        # NOTE: The chains given here refer to the `label_asym_id` field
        # of the `atom_site` category
        # However, by default `PDBxFile` uses the `auth_asym_id` as
        # chain ID
        matrices = []
        for id, op_expr, asym_id_expr in zip(
            assembly_gen_category["assembly_id"],
            assembly_gen_category["oper_expression"],
            assembly_gen_category["asym_id_list"],
        ):
            # Find the operation expressions for given assembly ID
            # We already asserted that the ID is actually present
            if id == assembly_id:
                operations = _parse_operation_expression(op_expr)
                affected_chain_ids = asym_id_expr.split(",")
                for i, operation in enumerate(operations):
                    rotations = []
                    translations = [] 
                    for op_step in operation:
                        rotation, translation = transformation_dict[op_step]
                        rotations.append(rotation)
                        translations.append(translation)
                    matrix = _chain_transformations(rotations, translations)
                    matrices.append((affected_chain_ids, matrix.tolist()))

        return matrices
    
    def get_assemblies(self):
        assembly_dict = {}
        for assembly_id in self.list_assemblies():
            assembly_dict[assembly_id] = self.get_transformations(assembly_id)
        
        return assembly_dict


def _chain_transformations(rotations, translations):
    """
    Get a total rotation/translation transformation by combining
    multiple rotation/translation transformations.
    This is done by intermediately combining rotation matrices and
    translation vectors into 4x4 matrices in the form
    
    |r11 r12 r13 t1|
    |r21 r22 r23 t2|
    |r31 r32 r33 t3|
    |0   0   0   1 |.
    """
    total_matrix = np.identity(4)
    for rotation, translation in zip(rotations, translations):
        matrix = np.zeros((4,4))
        matrix[:3, :3] = rotation
        matrix[:3, 3] = translation
        matrix[3, 3] = 1
        total_matrix = matrix @ total_matrix
    
    # return total_matrix[:3, :3], total_matrix[:3, 3]
    return matrix



def _get_transformations(struct_oper):
    """
    Get transformation operation in terms of rotation matrix and
    translation for each operation ID in ``pdbx_struct_oper_list``.
    """
    transformation_dict = {}
    for index, id in enumerate(struct_oper["id"]):
        rotation_matrix = np.array(
            [
                [
                    float(struct_oper[f"matrix[{i}][{j}]"][index])
                    for j in (1, 2, 3)
                ]
                for i in (1, 2, 3)
            ]
        )
        translation_vector = np.array(
            [float(struct_oper[f"vector[{i}]"][index]) for i in (1, 2, 3)]
        )
        transformation_dict[id] = (rotation_matrix, translation_vector)
    return transformation_dict


def _parse_operation_expression(expression):
    """
    Get successive operation steps (IDs) for the given
    ``oper_expression``.
    Form the cartesian product, if necessary.
    """
    # Split groups by parentheses:
    # use the opening parenthesis as delimiter
    # and just remove the closing parenthesis
    expressions_per_step = expression.replace(")", "").split("(")
    expressions_per_step = [e for e in expressions_per_step if len(e) > 0]
    # Important: Operations are applied from right to left
    expressions_per_step.reverse()

    operations = []
    for expr in expressions_per_step:
        if "-" in expr:
            if "," in expr:
                for gexpr in expr.split(","):
                    if "-" in gexpr:
                        first, last = gexpr.split("-")
                        operations.append(
                            [str(id) for id in range(int(first), int(last) + 1)]
                        )
                    else:
                        operations.append([gexpr])
            else:
                # Range of operation IDs, they must be integers
                first, last = expr.split("-")
                operations.append(
                    [str(id) for id in range(int(first), int(last) + 1)]
                )
        elif "," in expr:
            # List of operation IDs
            operations.append(expr.split(","))
        else:
            # Single operation ID
            operations.append([expr])

    # Cartesian product of operations
    return list(itertools.product(*operations))
