import numpy as np
import itertools

from . import AssemblyParser


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
        transformations = []
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
                    total_rotation, total_translation = _chain_transformations(
                        rotations, translations
                    )
                    transformations.append((
                        np.array(affected_chain_ids, dtype="U4").tolist(),
                        total_rotation.tolist(),
                        total_translation.tolist()
                    ))
        
        return transformations
    
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
    
    return total_matrix[:3, :3], total_matrix[:3, 3]



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