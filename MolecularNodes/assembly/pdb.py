import numpy as np
# import biotite.structure.io.mmtf as mmtf
from . import AssemblyParser


class PDBAssemblyParser(AssemblyParser):
    ### Implementation adapted from ``biotite.structure.io.pdb.file``

    def __init__(self, pdb_file):
        self._file = pdb_file
    

    def list_assemblies(self):
        return self._file.list_assemblies()
    

    def get_transformations(self, assembly_id):
        import biotite
        # Get lines containing transformations for assemblies
        remark_lines = self._file.get_remark(350)
        if remark_lines is None:
            raise biotite.InvalidFileError(
                "File does not contain assembly information (REMARK 350)"
            )
        # Get lines corresponding to selected assembly ID
        assembly_start_i = None
        assembly_stop_i = None
        for i, line in enumerate(remark_lines):
            if line.startswith("BIOMOLECULE"):
                current_assembly_id = line[12:].strip()
                if assembly_start_i is not None:
                    # Start was already found -> this is the next entry
                    # -> this is the stop
                    assembly_stop_i = i
                    break
                if current_assembly_id == assembly_id:
                    assembly_start_i = i
        # In case of the final assembly of the file,
        # the 'stop' is the end of REMARK 350 lines
        assembly_stop_i = len(remark_lines) if assembly_stop_i is None else i
        if assembly_start_i is None:
            raise KeyError(
                f"The assembly ID '{assembly_id}' is not found"
            )
        assembly_lines = remark_lines[assembly_start_i : assembly_stop_i]

        # Get transformations for a sets of chains
        transformations = []
        chain_set_start_indices = [
            i for i, line in enumerate(assembly_lines)
            if line.startswith("APPLY THE FOLLOWING TO CHAINS")
        ]
        # Add exclusive stop at end of records
        chain_set_start_indices.append(len(assembly_lines))
        assembly = None
        for i in range(len(chain_set_start_indices) - 1):
            start = chain_set_start_indices[i]
            stop = chain_set_start_indices[i+1]
            # Read affected chain IDs from the following line(s)
            affected_chain_ids = []
            transform_start = None
            for j, line in enumerate(assembly_lines[start : stop]):
                if line.startswith("APPLY THE FOLLOWING TO CHAINS:") or \
                   line.startswith("                   AND CHAINS:"):
                        affected_chain_ids += [
                            chain_id.strip() 
                            for chain_id in line[30:].split(",")
                        ]
                else:
                    # Chain specification has finished
                    # BIOMT lines start directly after chain specification
                    transform_start = start + j
                    break
            # Parse transformations from BIOMT lines
            if transform_start is None:
                raise biotite.InvalidFileError(
                    "No 'BIOMT' records found for chosen assembly"
                )
            rotations, translations = _parse_transformations(
                assembly_lines[transform_start : stop]
            )
            for rotation, translation in zip(rotations, translations):
                transformations.append((
                    np.array(affected_chain_ids, dtype="U4").tolist(),
                    rotation.tolist(),
                    translation.tolist()
                ))

        return transformations
    
    def get_assemblies(self):
        assembly_dict = {}
        for assembly_id in self.list_assemblies():
            assembly_dict[assembly_id] = self.get_transformations(assembly_id)
        
        return assembly_dict


def _parse_transformations(lines):
    """
    Parse the rotation and translation transformations from
    *REMARK* 290 or 350.
    Return as array of matrices and vectors respectively
    """
    import biotite
    # Each transformation requires 3 lines for the (x,y,z) components
    if len(lines) % 3 != 0:
        raise biotite.InvalidFileError("Invalid number of transformation vectors")
    n_transformations = len(lines) // 3

    rotations = np.zeros((n_transformations, 3, 3), dtype=float)
    translations = np.zeros((n_transformations, 3), dtype=float)

    transformation_i = 0
    component_i = 0
    for line in lines:
        # The first two elements (component and
        # transformation index) are not used
        transformations = [float(e) for e in line.split()[2:]]
        if len(transformations) != 4:
            raise biotite.InvalidFileError(
                "Invalid number of transformation vector elements"
            )
        rotations[transformation_i, component_i, :] = transformations[:3]
        translations[transformation_i, component_i] = transformations[3]

        component_i += 1
        if component_i == 3:
            # All (x,y,z) components were parsed
            # -> head to the next transformation 
            transformation_i += 1
            component_i = 0
    
    return rotations, translations