import itertools
import warnings

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import numpy as np
from biotite import InvalidFileError

from .assembly import AssemblyParser
from .molecule import Molecule


class OldCIF(Molecule):
    def __init__(self, file_path, extra_fields=None, sec_struct=True):
        super().__init__(file_path=file_path)
        self.array = self._get_structure(
            extra_fields=extra_fields, sec_struct=sec_struct
        )
        self.n_atoms = self.array.array_length()

    def _read(self, file_path):
        return pdbx.legacy.PDBxFile.read(file_path)

    def _get_structure(self, extra_fields: str = None, sec_struct=True, bonds=True):
        fields = ["b_factor", "charge", "occupancy", "atom_id"]
        if extra_fields:
            [fields.append(x) for x in extra_fields]

        # if the 'atom_site' doesn't exist then it will just be a small molecule
        # which can be extracted with the get_component()
        try:
            array = pdbx.get_structure(self.file, extra_fields=extra_fields)
            annotations = {
                "sec_struct": _get_secondary_structure,
                "entity_id": _get_entity_id,
            }
            for key, func in annotations.items():
                try:
                    array.set_annotation(key, func(array, self.file))
                except KeyError as e:
                    print(e)

        except InvalidFileError:
            array = pdbx.get_component(self.file)

        # pdbx files don't seem to have bond information defined, so connect them based
        # on their residue names
        if not array.bonds and bonds:
            array.bonds = struc.bonds.connect_via_residue_names(
                array, inter_residue=True
            )

        return array

    def _entity_ids(self):
        entities = self.file["entity"]
        if not entities:
            return None

        return entities.get("pdbx_description", None)

    def _assemblies(self):
        return CIFAssemblyParser(self.file).get_assemblies()


def _ss_label_to_int(label):
    if "HELX" in label:
        return 1
    elif "STRN" in label:
        return 2
    else:
        return 3


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

    # get the annotations for the struc_conf cetegory. Provides start and end
    # residues for the annotations. For most files this will only contain the
    # alpha helices, but will sometimes contain also other secondary structure
    # information such as in AlphaFold predictions

    conf = file.get_category("struct_conf")
    if not conf:
        raise KeyError
    starts = conf["beg_auth_seq_id"].astype(int)
    ends = conf["end_auth_seq_id"].astype(int)
    chains = conf["end_auth_asym_id"].astype(str)
    id_label = conf["id"].astype(str)

    # most files will have a separate category for the beta sheets
    # this can just be appended to the other start / end / id and be processed
    # as normal
    sheet = file.get_category("struct_sheet_range")
    if sheet:
        starts = np.append(starts, sheet["beg_auth_seq_id"].astype(int))
        ends = np.append(ends, sheet["end_auth_seq_id"].astype(int))
        chains = np.append(chains, sheet["end_auth_asym_id"].astype(str))
        id_label = np.append(id_label, np.repeat("STRN", len(sheet["id"])))

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


def _get_entity_id(array, file):
    entities = file.get_category("entity_poly")
    if not entities:
        raise KeyError
    chain_ids = entities["pdbx_strand_id"]

    # the chain_ids are an array of individual items np.array(['A,B', 'C', 'D,E,F'])
    # which need to be categorised as [1, 1, 2, 3, 3, 3] for their belonging to individual
    # entities

    chains = []
    idx = []
    for i, chain_str in enumerate(chain_ids):
        for chain in chain_str.split(","):
            chains.append(chain)
            idx.append(i)

    entity_lookup = dict(zip(chains, idx))
    chain_id_int = np.array(
        [entity_lookup.get(chain, -1) for chain in array.chain_id], int
    )
    return chain_id_int


class CIFAssemblyParser(AssemblyParser):
    # Implementation adapted from ``biotite.structure.io.pdbx.convert``

    def __init__(self, file_cif):
        self._file = file_cif

    def list_assemblies(self):
        return list(pdbx.list_assemblies(self._file).keys())

    def get_transformations(self, assembly_id):
        assembly_gen_category = self._file["pdbx_struct_assembly_gen"]

        struct_oper_category = self._file["pdbx_struct_oper_list"]

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
    for index, id in enumerate(struct_oper["id"]):
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
