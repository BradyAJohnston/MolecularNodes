import itertools
from io import BytesIO
from pathlib import Path
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import numpy as np
from biotite import InvalidFileError
from .reader import ReaderBase


class PDBXReader(ReaderBase):
    def __init__(self, file_path: str | Path | BytesIO):
        self._extra_annotations = {
            "sec_struct": self._get_secondary_structure,
            "entity_id": self._get_entity_id,
        }
        self._extra_fields = ["b_factor", "occupancy", "atom_id"]
        super().__init__(file_path)

    def read(self, file_path):
        match file_path.suffix:
            case ".cif":
                return pdbx.CIFFile.read(file_path)
            case ".bcif":
                return pdbx.BinaryCIFFile.read(file_path)
            case _:
                raise NotImplementedError(
                    f"File type {file_path.suffix} not supported."
                )

    def get_structure(
        self, model: int | None = None
    ) -> struc.AtomArray | struc.AtomArrayStack:
        try:
            array = pdbx.get_structure(
                self.file, model=model, extra_fields=self._extra_fields
            )
        except InvalidFileError:
            array = pdbx.get_component(self.file)

        if not array.bonds:
            array.bonds = struc.bonds.connect_via_residue_names(
                array, inter_residue=True
            )

        return array  # type: ignore

    def _assemblies(self):
        return CIFAssemblyParser(self.file).get_assemblies()

    def entity_ids(self):
        try:
            return (
                self.file.block.get("entity")
                .get("pdbx_description")
                .as_array()
                .tolist()
            )
        except AttributeError:
            return None

    @staticmethod
    def _extract_matrices(category):
        matrix_columns = [
            "matrix[1][1]",
            "matrix[1][2]",
            "matrix[1][3]",
            "vector[1]",
            "matrix[2][1]",
            "matrix[2][2]",
            "matrix[2][3]",
            "vector[2]",
            "matrix[3][1]",
            "matrix[3][2]",
            "matrix[3][3]",
            "vector[3]",
        ]

        columns = [category[name].as_array().astype(float) for name in matrix_columns]
        matrices = np.empty((len(columns[0]), 4, 4), float)

        col_mask = np.tile((0, 1, 2, 3), 3)
        row_mask = np.repeat((0, 1, 2), 4)
        for column, coli, rowi in zip(columns, col_mask, row_mask):
            matrices[:, rowi, coli] = column

        return matrices

    @staticmethod
    def _get_entity_id(array, file):
        chain_ids = file.block["entity_poly"]["pdbx_strand_id"].as_array(str)

        # the chain_ids are an array of individual items np.array(['A,B', 'C', 'D,E,F'])
        # which need to be categorised as [1, 1, 2, 3, 3, 3] for their belonging to individual
        # entities

        chains = []
        idx = []
        for i, chain_str in enumerate(chain_ids):
            for chain in chain_str.split(","):
                chains.append(chain)
                idx.append(i)

        # this is how we map the chain_ids and res_names of our entities to their integer
        # representations
        entity_lookup = dict(zip(chains, idx))

        # for the hetero atoms, we need to add a new entity_id into the lookup so that
        # they can be assigned an entity ID
        unique_res_het = np.unique(array.res_name[array.hetero])
        for het in unique_res_het:
            if het not in entity_lookup:
                entity_lookup[het] = max(entity_lookup.values()) + 1

        entity_id_int = np.zeros(len(array.res_name), int)
        for i, res_name in enumerate(array.res_name):
            entity_id_int[i] = entity_lookup.get(
                res_name, entity_lookup[array.chain_id[i]]
            )

        return entity_id_int

    @staticmethod
    def _get_secondary_structure(array, file):
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

        # get the annotations for the struc_conf category. Provides start and end
        # residues for the annotations. For most files this will only contain the
        # alpha helices, but will sometimes contain also other secondary structure
        # information such as in AlphaFold predictions

        conf = file.block.get("struct_conf")
        if conf:
            starts = conf["beg_auth_seq_id"].as_array().astype(int)
            ends = conf["end_auth_seq_id"].as_array().astype(int)
            chains = conf["end_auth_asym_id"].as_array().astype(str)
            id_label = conf["id"].as_array().astype(str)
        else:
            starts = np.empty(0, dtype=int)
            ends = np.empty(0, dtype=int)
            chains = np.empty(0, dtype=str)
            id_label = np.empty(0, dtype=int)

        # most files will have a separate category for the beta sheets
        # this can just be appended to the other start / end / id and be processed
        # as normalquit
        sheet = file.block.get("struct_sheet_range")
        if sheet:
            starts = np.append(starts, sheet["beg_auth_seq_id"].as_array().astype(int))
            ends = np.append(ends, sheet["end_auth_seq_id"].as_array().astype(int))
            chains = np.append(chains, sheet["end_auth_asym_id"].as_array().astype(str))
            id_label = np.append(id_label, np.repeat("STRN", len(sheet["id"])))

        if not conf and not sheet:
            raise KeyError

        # convert the string labels to integer representations of the SS
        # AH: 1, BS: 2, LOOP: 3

        id_int = np.array([_ss_label_to_int(label) for label in id_label], int)

        # create a lookup dictionary that enables lookup of secondary structure
        # based on the chain_id and res_id values

        lookup = dict()
        for chain in np.unique(chains):
            arrays = []
            mask = chain == chains
            start_sub = starts[mask]
            end_sub = ends[mask]
            id_sub = id_int[mask]

            for start, end, id in zip(start_sub, end_sub, id_sub):
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

    for group in oper.strip(")").split("("):
        if "," in group:
            for i in group.split(","):
                op_ids.append()

    for group in oper.split(","):
        if "-" not in group:
            op_ids.append(str(group))
            continue

        start, stop = [int(x) for x in group.strip("()").split("-")]
        for i in range(start, stop + 1):
            op_ids.append(str(i))

    return op_ids


def _ss_label_to_int(label):
    if "HELX" in label:
        return 1
    elif "STRN" in label:
        return 2
    else:
        return 3


class CIFAssemblyParser:
    # Implementation adapted from ``biotite.structure.io.pdbx.convert``

    def __init__(self, file_cif):
        self._file = file_cif

    def list_assemblies(self):
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
        pdb_model_num = -1
        for id, op_expr, asym_id_expr in zip(
            assembly_gen_category["assembly_id"].as_array(str),
            assembly_gen_category["oper_expression"].as_array(str),
            assembly_gen_category["asym_id_list"].as_array(str),
        ):
            pdb_model_num += 1
            # Find the operation expressions for given assembly ID
            # We already asserted that the ID is actually present
            if id != assembly_id:
                continue

            operations = _parse_operation_expression(op_expr)

            affected_chain_ids = asym_id_expr.split(",")

            for i, operation in enumerate(operations):
                # for op_step in operation:
                matrices.append(
                    {
                        "chain_ids": affected_chain_ids,
                        "matrix": transformation_dict[operation[0]].tolist(),
                        "pdb_model_num": pdb_model_num,
                    }
                )

        return matrices

    def get_assemblies(self):
        assembly_dict = {}
        for assembly_id in self.list_assemblies():
            assembly_dict[assembly_id] = self.get_transformations(assembly_id)

        return assembly_dict


def _extract_matrices(category, scale=True):
    matrix_columns = [
        "matrix[1][1]",
        "matrix[1][2]",
        "matrix[1][3]",
        "vector[1]",
        "matrix[2][1]",
        "matrix[2][2]",
        "matrix[2][3]",
        "vector[2]",
        "matrix[3][1]",
        "matrix[3][2]",
        "matrix[3][3]",
        "vector[3]",
    ]

    columns = [category[name].as_array().astype(float) for name in matrix_columns]
    n = 4 if scale else 3
    matrices = np.empty((len(columns[0]), n, 4), float)

    col_mask = np.tile((0, 1, 2, 3), 3)
    row_mask = np.repeat((0, 1, 2), 4)
    for column, coli, rowi in zip(columns, col_mask, row_mask):
        matrices[:, rowi, coli] = column

    return dict(zip(category["id"].as_array(str), matrices))


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
                [float(struct_oper[f"matrix[{i}][{j}]"][index]) for j in (1, 2, 3)]
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
                operations.append([str(id) for id in range(int(first), int(last) + 1)])
        elif "," in expr:
            # List of operation IDs
            operations.append(expr.split(","))
        else:
            # Single operation ID
            operations.append([expr])

    # Cartesian product of operations
    return list(itertools.product(*operations))
