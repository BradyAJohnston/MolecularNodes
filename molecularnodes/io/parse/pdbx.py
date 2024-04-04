import numpy as np
import warnings
import itertools

from .molecule import Molecule


class PDBX(Molecule):
    def __init__(self, file_path):
        self.file_path = file_path
        self.file = self._read(file_path)

    @property
    def entity_ids(self):
        return self.file.block.get('entity').get('pdbx_description').as_array().tolist()

    def _get_entity_id(self, array, file):
        chain_ids = file.block['entity_poly']['pdbx_strand_id'].as_array()

        # the chain_ids are an array of individual items np.array(['A,B', 'C', 'D,E,F'])
        # which need to be categorised as [1, 1, 2, 3, 3, 3] for their belonging to individual
        # entities

        chains = []
        idx = []
        for i, chain_str in enumerate(chain_ids):
            for chain in chain_str.split(','):
                chains.append(chain)
                idx.append(i)

        entity_lookup = dict(zip(chains, idx))
        chain_id_int = np.array([entity_lookup.get(chain, -1)
                                for chain in array.chain_id], int)
        return chain_id_int

    def get_structure(self, extra_fields=['b_factor', 'occupancy', 'atom_id'], bonds=True):
        import biotite.structure.io.pdbx as pdbx
        import biotite.structure as struc

        array = pdbx.get_structure(self.file, extra_fields=extra_fields)
        try:
            array.set_annotation(
                'sec_struct', self._get_secondary_structure(
                    array=array, file=self.file)
            )
        except KeyError:
            warnings.warn('No secondary structure information.')
        try:
            array.set_annotation(
                'entity_id', self._get_entity_id(array, self.file)
            )
        except KeyError:
            warnings.warn('No entity ID information')

        if not array.bonds and bonds:
            array.bonds = struc.bonds.connect_via_residue_names(
                array, inter_residue=True)

        return array

    def _assemblies(self):
        return CIFAssemblyParser(self.file).get_assemblies()

        # # in the cif / BCIF file 3x4 transformation matrices are stored in individual
        # # columns, this extracts them and returns them with additional row for scaling,
        # # meaning an (n, 4, 4) array is returned, where n is the number of transformations
        # # and each is a 4x4 transformaiton matrix
        # cat_matrix = self.file.block['pdbx_struct_oper_list']
        # matrices = self._extract_matrices(cat_matrix)

        # # sometimes there will be missing opers / matrices. For example in the
        # # 'square.bcif' file, the matrix IDs go all the way up to 18024, but only
        # # 18023 matrices are defined. That is becuase matrix 12 is never referenced, so
        # # isn't included in teh file. To get around this we have to just get the specific
        # # IDs that are defined for the matrices and use that to lookup the correct index
        # # in the matrices array.
        # mat_ids = cat_matrix.get('id').as_array(int)
        # mat_lookup = dict(zip(mat_ids, range(len(mat_ids))))

        # category = self.file.block['pdbx_struct_assembly_gen']
        # ids = category['assembly_id'].as_array(int)
        # opers = category['oper_expression'].as_array(str)
        # asyms = category['asym_id_list'].as_array()

        # # constructs a dictionary of
        # # {
        # #   '1': ((['A', 'B', C'], [4x4 matrix]), (['A', 'B'], [4x4 matrix])),
        # #   '2': ((['A', 'B', C'], [4x4 matrix]))
        # # }
        # # where each entry in the dictionary is a biological assembly, and each dictionary
        # # value contains a list of tranasformations which need to be applied. Each entry in
        # # the list of transformations is
        # # ([chains to be affected], [4x4 transformation matrix])
        # assembly_dic = {}
        # for idx, oper, asym in zip(ids, opers, asyms):
        #     trans = list()
        #     asym = asym.split(',')
        #     for op in _parse_opers(oper):
        #         i = int(op)
        #         trans.append((asym, matrices[mat_lookup[i]].tolist()))
        #     assembly_dic[str(idx)] = trans

        # return assembly_dic

    def _extract_matrices(self, category):
        matrix_columns = [
            'matrix[1][1]',
            'matrix[1][2]',
            'matrix[1][3]',
            'vector[1]',
            'matrix[2][1]',
            'matrix[2][2]',
            'matrix[2][3]',
            'vector[2]',
            'matrix[3][1]',
            'matrix[3][2]',
            'matrix[3][3]',
            'vector[3]'
        ]

        columns = [
            category[name].as_array().astype(float) for
            name in matrix_columns
        ]
        matrices = np.empty((len(columns[0]), 4, 4), float)

        col_mask = np.tile((0, 1, 2, 3), 3)
        row_mask = np.repeat((0, 1, 2), 4)
        for column, coli, rowi in zip(columns, col_mask, row_mask):
            matrices[:, rowi, coli] = column

        return matrices

    def _get_secondary_structure(self, file, array):
        """
        Get secondary structure information for the array from the file.

        Parameters
        ----------
        array : numpy array
            The array for which secondary structure information is to be retrieved.
        file : object
            The file object containing the secondary structure information.

        Returns
        -------
        numpy array
            A numpy array of secondary structure information, where each element is either 0, 1, 2, or 3.
            - 0: Not a peptide
            - 1: Alpha helix
            - 2: Beta sheet
            - 3: Loop

        Raises
        ------
        KeyError
            If the 'struct_conf' category is not found in the file.
        """
        import biotite.structure as struc

        # get the annotations for the struc_conf cetegory. Provides start and end
        # residues for the annotations. For most files this will only contain the
        # alpha helices, but will sometimes contain also other secondary structure
        # information such as in AlphaFold predictions

        conf = file.block.get('struct_conf')
        if not conf:
            raise KeyError
        starts = conf['beg_auth_seq_id'].as_array().astype(int)
        ends = conf['end_auth_seq_id'].as_array().astype(int)
        chains = conf['end_auth_asym_id'].as_array().astype(str)
        id_label = conf['id'].as_array().astype(str)

        # most files will have a separate category for the beta sheets
        # this can just be appended to the other start / end / id and be processed
        # as normalquit
        sheet = file.block.get('struct_sheet_range')
        if sheet:
            starts = np.append(
                starts, sheet['beg_auth_seq_id'].as_array().astype(int))
            ends = np.append(
                ends, sheet['end_auth_seq_id'].as_array().astype(int))
            chains = np.append(
                chains, sheet['end_auth_asym_id'].as_array().astype(str))
            id_label = np.append(id_label, np.repeat('STRN', len(sheet['id'])))

        # convert the string labels to integer representations of the SS
        # AH: 1, BS: 2, LOOP: 3

        id_int = np.array([_ss_label_to_int(label) for label in id_label], int)

        # create a lookup dictionary that enables lookup of secondary structure
        # based on the chain_id and res_id values

        lookup = dict()
        for chain in np.unique(chains):
            arrays = []
            mask = (chain == chains)
            start_sub = starts[mask]
            end_sub = ends[mask]
            id_sub = id_int[mask]

            for (start, end, id) in zip(start_sub, end_sub, id_sub):
                idx = np.arange(start, end + 1, dtype=int)
                arr = np.zeros((len(idx), 2), dtype=int)
                arr[:, 0] = idx
                arr[:, 1] = 3
                arr[:, 1] = id
                arrays.append(arr)

            lookup[chain] = dict(np.vstack(arrays).tolist())

        # use the lookup dictionary to get the SS annotation based on the chain_id and res_id
        secondary_structure = np.zeros(len(array.chain_id), int)
        for i, (chain, res) in enumerate(zip(array.chain_id, array.res_id)):
            try:
                secondary_structure[i] = lookup[chain].get(res, 3)
            except KeyError:
                secondary_structure[i] = 0

        # assign SS to 0 where not peptide
        secondary_structure[~struc.filter_amino_acids(array)] = 0
        return secondary_structure


def _parse_opers(oper):
    # we want the example '1,3,(5-8)' to expand to (1, 3, 5, 6, 7, 8).
    op_ids = list()

    for group in oper.strip(')').split('('):
        if "," in group:
            for i in group.split(','):
                op_ids.append()

    for group in oper.split(","):
        if "-" not in group:
            op_ids.append(str(group))
            continue

        start, stop = [int(x) for x in group.strip("()").split('-')]
        for i in range(start, stop + 1):
            op_ids.append(str(i))

    return op_ids


def _ss_label_to_int(label):
    if 'HELX' in label:
        return 1
    elif 'STRN' in label:
        return 2
    else:
        return 3


class CIF(PDBX):
    def __init__(self, file_path):
        super().__init__(file_path)
        self.file_path = file_path
        # self.file = self.read(file_path)
        self.array = self.get_structure()

    def _read(self, file_path):
        import biotite.structure.io.pdbx as pdbx
        return pdbx.CIFFile.read(file_path)


class BCIF(PDBX):
    def __init__(self, file_path):
        super().__init__(file_path)
        self.file_path = file_path
        # self.file = self.read(file_path)
        self.array = self.get_structure()

    def _read(self, file_path):
        import biotite.structure.io.pdbx as pdbx
        return pdbx.BinaryCIFFile.read(file_path)


class CIFAssemblyParser:
    # Implementation adapted from ``biotite.structure.io.pdbx.convert``

    def __init__(self, file_cif):
        self._file = file_cif

    def list_assemblies(self):
        import biotite.structure.io.pdbx as pdbx
        return list(pdbx.list_assemblies(self._file).keys())

    def get_transformations(self, assembly_id):
        assembly_gen_category = self._file.block["pdbx_struct_assembly_gen"]

        struct_oper_category = self._file.block["pdbx_struct_oper_list"]

        if assembly_id not in assembly_gen_category["assembly_id"].as_array(str):
            raise KeyError(f"File has no Assembly ID '{assembly_id}'")

        # Extract all possible transformations indexed by operation ID
        # transformation_dict = _get_transformations(struct_oper_category)
        transformation_dict = _extract_matrices(struct_oper_category)

        # Get necessary transformations and the affected chain IDs
        # NOTE: The chains given here refer to the `label_asym_id` field
        # of the `atom_site` category
        # However, by default `PDBxFile` uses the `auth_asym_id` as
        # chain ID
        matrices = []
        for id, op_expr, asym_id_expr in zip(
            assembly_gen_category["assembly_id"].as_array(str),
            assembly_gen_category["oper_expression"].as_array(str),
            assembly_gen_category["asym_id_list"].as_array(str),
        ):
            # Find the operation expressions for given assembly ID
            # We already asserted that the ID is actually present
            if id == assembly_id:
                operations = _parse_operation_expression(op_expr)
                affected_chain_ids = asym_id_expr.split(",")
                for i, operation in enumerate(operations):
                    for op_step in operation:
                        matrix = transformation_dict[op_step]
                    matrices.append((affected_chain_ids, matrix.tolist()))

        return matrices

    def get_assemblies(self):
        assembly_dict = {}
        for assembly_id in self.list_assemblies():
            assembly_dict[assembly_id] = self.get_transformations(assembly_id)

        return assembly_dict


def _extract_matrices(category, scale=True):
    matrix_columns = [
        'matrix[1][1]',
        'matrix[1][2]',
        'matrix[1][3]',
        'vector[1]',
        'matrix[2][1]',
        'matrix[2][2]',
        'matrix[2][3]',
        'vector[2]',
        'matrix[3][1]',
        'matrix[3][2]',
        'matrix[3][3]',
        'vector[3]'
    ]

    columns = [
        category[name].as_array().astype(float) for
        name in matrix_columns
    ]
    n = 4 if scale else 3
    matrices = np.empty((len(columns[0]), n, 4), float)

    col_mask = np.tile((0, 1, 2, 3), 3)
    row_mask = np.repeat((0, 1, 2), 4)
    for column, coli, rowi in zip(columns, col_mask, row_mask):
        matrices[:, rowi, coli] = column

    return dict(zip(category['id'].as_array(str), matrices))


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
        matrix = np.zeros((4, 4))
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
    for index, id in enumerate(struct_oper["id"].as_array()):
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
                            [str(id)
                             for id in range(int(first), int(last) + 1)]
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
