import numpy as np
from biotite.structure import AtomArray, connect_via_distances
from ...assets import data
from ..density.cube import read_cube
from .reader import ReaderBase


class CubeReader(ReaderBase):
    """Read the atoms of a Gaussian cube file, guessing bonds from distances."""

    _extra_annotations = {}

    def read(self, file_path):
        return read_cube(file_path)

    def get_structure(self):
        cube = self.file
        elements = [
            data.elements_by_atomic_number.get(int(z), {}).get("element_symbol", "X")
            for z in cube.atomic_numbers
        ]
        array = AtomArray(len(elements))
        array.coord = cube.positions
        array.element = np.char.upper(elements)
        array.atom_name = np.char.upper(elements)
        array.res_name[:] = "UNL"
        array.res_id[:] = 1
        array.chain_id[:] = "A"
        array.hetero[:] = True
        array.bonds = connect_via_distances(array, inter_residue=False)
        return array

    def _assemblies(self):
        return None
