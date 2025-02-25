import numpy as np
import biotite.structure as struc
from biotite import InvalidFileError
from biotite.structure import (
    BadStructureError,
    annotate_sse,
    spread_residue_wise,
    connect_via_residue_names,
    AtomArray
)

import uuid
from .base import Molecule
from .pdb import _comp_secondary_structure
from ..base import MolecularEntity, EntityType



class Array(Molecule):
    def __init__(self, array: AtomArray, name="AtomArray"):
        self.file_path = None
        self.file = "ARRAY_LOADED_DIRECTLY"
        self._frames_collection = None
        self._entity_type = EntityType.MOLECULE
        self._assemblies = lambda: None
        self._uuid = str(uuid.uuid4())
        self.array = self._validate_structure(array)
        self.n_atoms = self.array.array_length()    # note: I am not sure if this is consistent
        self.create_object(name=name)               # Create the Blender object

    def read(self, file_path):
        pass

    def _validate_structure(self, array: AtomArray):
        print(array.get_annotation_categories())
        # TODO: implement entity ID, sec_struct for PDB files
        extra_fields = ["b_factor", "occupancy", "charge", "atom_id"]

        # https://github.com/biotite-dev/biotite/blob/main/src/biotite/structure/io/pdb/file.py#L331
        # for field in extra_fields:
        #     if field == "atom_id":
        #         # Copy is necessary to avoid double masking in
        #         # later altloc ID filtering
        #         # array.set_annotation("atom_id", atom_id.copy())
        #         pass
        #     elif field == "charge":
        #         charge = np.array(charge_raw)
        #         array.set_annotation(
        #             "charge", np.where(charge == "  ", "0", charge).astype(int)
        #         )
        #     elif field == "occupancy":
        #         array.set_annotation("occupancy", occupancy)
        #     elif field == "b_factor":
        #         array.set_annotation("b_factor", b_factor)
        #     else:
        #         raise ValueError(f"Unknown extra field: {field}")

        sec_struct = _comp_secondary_structure(array)
        array.set_annotation("sec_struct", sec_struct)
        return array

    def _assemblies(self):
        return None
