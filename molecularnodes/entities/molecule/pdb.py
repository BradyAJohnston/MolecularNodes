import biotite.structure as struc
import numpy as np
from biotite import InvalidFileError
from biotite.structure import (
    BadStructureError,
    annotate_sse,
    connect_via_residue_names,
    spread_residue_wise,
)
from biotite.structure.io import pdb
from .assembly import AssemblyParser
from .reader import ReaderBase


class PDBReader(ReaderBase):
    extra_fields = ["b_factor", "occupancy", "charge", "atom_id"]
    _extra_annotations = {}

    @staticmethod
    def set_extra_annotations(annotations, array, file):
        return array

    def read(self, file_path):
        return pdb.PDBFile.read(file_path)

    def get_structure(self, model: int | None = None):
        # a bit dirty, but we first try and get the bond information from the file
        # if that fails, then we extract without the bonds and try to create bonds based
        # on residue / atom names.
        try:
            array = pdb.get_structure(
                pdb_file=self.file,
                model=model,
                extra_fields=self.extra_fields,
                include_bonds=True,
            )
        except AttributeError as e:
            print(
                f"Unable to get bond information: {e}\nAttempting `connect_via_residue_names()`"
            )
            array = pdb.get_structure(
                pdb_file=self.file,
                model=model,
                extra_fields=self.extra_fields,
                include_bonds=False,
            )
            try:
                array.bonds = connect_via_residue_names(array)
            except AttributeError as e:
                print("Not able to find bonds via residue: {e}")

        try:
            sec_struct = _get_sec_struct(self.file, array)
        except BadStructureError:
            sec_struct = _comp_secondary_structure(array[0])

        array.set_annotation("sec_struct", sec_struct)

        if array is None:
            raise InvalidFileError("Unable to read PDB file.")

        return array

    def _assemblies(self):
        return PDBAssemblyParser(self.file).get_assemblies()


def _get_sec_struct(file, array):
    lines = np.array(file.lines)
    lines_helix = lines[np.char.startswith(lines, "HELIX")]
    lines_sheet = lines[np.char.startswith(lines, "SHEET")]
    if len(lines_helix) == 0 and len(lines_sheet) == 0:
        raise struc.BadStructureError("No secondary structure information detected.")

    sec_struct = np.zeros(array.array_length(), int)

    helix_values = (22, 25, 34, 37, 20)
    sheet_values = (23, 26, 34, 37, 22)

    values = ((lines_helix, 1, helix_values), (lines_sheet, 2, sheet_values))

    def _get_mask(line, start1, end1, start2, end2, chainid):
        """
        Takes a line and makes a mask for the array for atoms which are defined inside the
        range of values defined by the start and end values.
        """
        # bump the starting values down by one for indexing into
        # pythong strings
        start1 -= 1
        start2 -= 1
        chainid -= 1

        # get the values as intergers
        start_num = int(line[start1:end1])
        end_num = int(line[start2:end2])
        chain_id = line[chainid].strip()

        # create a mask for the array based on these values
        mask = np.logical_and(
            np.logical_and(array.chain_id == chain_id, array.res_id >= start_num),
            array.res_id <= end_num,
        )

        return mask

    # assigns secondary structure based on lines and stores it in the sec_struct
    for lines, idx, value_list in values:
        for line in lines:
            sec_struct[_get_mask(line, *value_list)] = idx

    # assign remaining AA atoms to 3 (loop), while all other remaining
    # atoms will be 0 (not relevant)
    mask = np.logical_and(sec_struct == 0, struc.filter_canonical_amino_acids(array))

    sec_struct[mask] = 3

    return sec_struct


def _comp_secondary_structure(array):
    """Use dihedrals to compute the secondary structure of proteins

    Through biotite built-in method derivated from P-SEA algorithm (Labesse 1997)
    Returns an array with secondary structure for each atoms where:
    - 0 = '' = non-protein or not assigned by biotite annotate_sse
    - 1 = a = alpha helix
    - 2 = b = beta sheet
    - 3 = c = coil

    Inspired from https://www.biotite-python.org/examples/gallery/structure/transketolase_sse.html
    """
    # TODO Port [PyDSSP](https://github.com/ShintaroMinami/PyDSSP)

    conv_sse_char_int = {"a": 1, "b": 2, "c": 3, "": 0}

    char_sse = annotate_sse(array)
    int_sse = np.array([conv_sse_char_int[char] for char in char_sse], dtype=int)
    atom_sse = spread_residue_wise(array, int_sse)

    return atom_sse


class PDBAssemblyParser(AssemblyParser):
    # Implementation adapted from ``biotite.structure.io.pdb.file``

    def __init__(self, pdb_file):
        self._file = pdb_file

    def list_assemblies(self):
        return self._file.list_assemblies()

    def get_transformations(self, assembly_id):
        # Get lines containing transformations for assemblies
        remark_lines = self._file.get_remark(350)
        if remark_lines is None:
            raise InvalidFileError(
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
            raise KeyError(f"The assembly ID '{assembly_id}' is not found")
        assembly_lines = remark_lines[assembly_start_i:assembly_stop_i]

        # Get transformations for a sets of chains
        transformations = []
        chain_set_start_indices = [
            i
            for i, line in enumerate(assembly_lines)
            if line.startswith("APPLY THE FOLLOWING TO CHAINS")
        ]
        # Add exclusive stop at end of records
        chain_set_start_indices.append(len(assembly_lines))
        for i in range(len(chain_set_start_indices) - 1):
            start = chain_set_start_indices[i]
            stop = chain_set_start_indices[i + 1]
            # Read affected chain IDs from the following line(s)
            affected_chain_ids = []
            transform_start = None
            for j, line in enumerate(assembly_lines[start:stop]):
                if line.startswith("APPLY THE FOLLOWING TO CHAINS:") or line.startswith(
                    "                   AND CHAINS:"
                ):
                    affected_chain_ids += [
                        chain_id.strip() for chain_id in line[30:].split(",")
                    ]
                else:
                    # Chain specification has finished
                    # BIOMT lines start directly after chain specification
                    transform_start = start + j
                    break
            # Parse transformations from BIOMT lines
            if transform_start is None:
                raise InvalidFileError("No 'BIOMT' records found for chosen assembly")

            matrices = _parse_transformations(assembly_lines[transform_start:stop])

            for matrix in matrices:
                transformations.append((affected_chain_ids, matrix.tolist()))

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

    # Each transformation requires 3 lines for the (x,y,z) components
    if len(lines) % 3 != 0:
        raise InvalidFileError("Invalid number of transformation vectors")
    n_transformations = len(lines) // 3

    matrices = np.tile(np.identity(4), (n_transformations, 1, 1))

    transformation_i = 0
    component_i = 0
    for line in lines:
        # The first two elements (component and
        # transformation index) are not used
        transformations = [float(e) for e in line.split()[2:]]

        if len(transformations) != 4:
            raise InvalidFileError("Invalid number of transformation vector elements")
        matrices[transformation_i, component_i, :] = transformations

        component_i += 1
        if component_i == 3:
            # All (x,y,z) components were parsed
            # -> head to the next transformation
            transformation_i += 1
            component_i = 0

    return matrices
