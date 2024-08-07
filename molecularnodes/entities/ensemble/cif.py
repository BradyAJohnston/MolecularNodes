from pathlib import Path
import numpy as np
from mathutils import Matrix
from typing import Any, Dict, List, Optional, TypedDict, Union
from biotite.structure import AtomArray
import biotite.structure.io.pdbx as pdbx


class CIF:
    def __init__(self, file_path):
        # super().__init__()
        self.file_path = file_path
        self.file = self.read()
        self.entities = {}
        categories = self.file.block
        # check if a petworld CellPack model or not
        self.is_petworld = False
        if "PDB_model_num" in categories["pdbx_struct_assembly_gen"]:
            print("PetWorld!")
            self.is_petworld = True
        entity = {}
        entityids = []
        pdbx_description = []
        if self.is_petworld:
            entity = categories['pdbx_model']
            entityids = [str(i+1) for i in range(len(entity['name']))]
            pdbx_description = entity['name'].as_array()
        else:
            entity = categories['entity']
            entityids = entity['id'].as_array()
            pdbx_description = entity['pdbx_description'].as_array()
        for i in range(len(entityids)):
            self.entities[entityids[i]] = pdbx_description[i]

        self.array = _atom_array_from_cif(categories)
        self._transforms_data = _get_ops_from_cif(categories)
        self.n_models = 1
        self.n_atoms = self.array.shape
        if self.is_petworld:
            self.array.asym_id = self.array.chain_id
            # self.array.chain_id = self.array.asym_id
        self.chain_ids = self._chain_ids()

    def read(self):
        suffix = Path(self.file_path).suffix
        print('reading file', self.file_path)
        if suffix in (".bin", ".bcif"):
            return pdbx.BinaryCIFFile.read(self.file_path)
        elif suffix == ".cif":
            return pdbx.CIFFile.read(self.file_path)
        # with open(self.file_path, "rb") as data:
        #     open_bcif = loads(data.read())
        #
        # return open_bcif

    def assemblies(self, as_array=True):
        return self._transforms_data

    def _chain_ids(self, as_int=False):
        if as_int:
            return np.unique(self.array.chain_id, return_inverse=True)[1]
        return np.unique(self.array.chain_id)


def _atom_array_from_cif(categories):
    # categories = open_bcif.data_blocks[0]

    # check if a petworld CellPack model or not
    is_petworld = False
    if "PDB_model_num" in categories["pdbx_struct_assembly_gen"]:
        print("PetWorld!")
        is_petworld = True

    atom_site = categories["atom_site"]
    n_atoms = atom_site.row_count

    # Initialise the atom array that will contain all of the data for the atoms
    # in the bcif file. TODO support multi-model bcif files
    # we first pull out the coordinates
    # as they are from 3 different fields, but all
    # other fields should be single self-contained fields
    mol = AtomArray(n_atoms)
    coord_field_names = [f"Cartn_{axis}" for axis in "xyz"]
    mol.coord = np.hstack(
        list(
            [
                atom_site[column].as_array().reshape((n_atoms, 1))
                for column in coord_field_names
            ]
        )
    )

    # the list of current
    atom_site_lookup = {
        # have to make sure the chain_id
        # ends up being the same as the space operator
        "label_asym_id": "chain_id",
        "label_atom_id": "atom_name",
        "label_comp_id": "res_name",
        "type_symbol": "element",
        "label_seq_id": "res_id",
        "B_iso_or_equiv": "b_factor",
        "label_entity_id": "entity_id",
        "pdbx_PDB_model_num": "model_id",
        "pdbx_formal_charge": "charge",
        "occupancy": "occupany",
        "id": "atom_id",
    }

    if is_petworld:
        print("PetWorld!")
        # annotations[0][1] = 'pdbx_PDB_model_num'
        atom_site_lookup.pop("label_asym_id")
        atom_site_lookup["pdbx_PDB_model_num"] = "chain_id"
        atom_site_lookup.pop("label_entity_id")

    # for name in atom_site.field_names:
    for name, column in atom_site.items():
        # the coordinates have already been extracted
        # so we can skip over those field names
        if name in coord_field_names:
            continue
        # numpy does a pretty good job of guessing
        # the data types from the fields
        data = atom_site[name].as_array()
        if name == "label_asym_id":
            # print("set annoatation ", name)
            # print(data)
            mol.asym_id = data
        # if a specific name for an annotation is
        # already specified earlier, we can
        # use that to ensure consitency. All other
        # fields are also still added as we
        # may as well do so, in case we want any extra data
        annotation_name = atom_site_lookup.get(name)
        if not annotation_name:
            annotation_name = name
        # TODO this could be expanded to capture
        # fields that are entirely '' and drop them
        # or fill them with 0s
        if annotation_name == "res_id" and (data[0] == "" or data[0] == "."):
            data = np.array([0 if (x == "" or x == ".") else x for x in data])
        mol.set_annotation(annotation_name, data)
        if name == "pdbx_PDB_model_num" and is_petworld:
            mol.set_annotation('entity_id', data)

    return mol


def rotation_from_matrix(matrix):
    rotation_matrix = np.identity(4, dtype=float)
    rotation_matrix[:3, :3] = matrix
    translation, rotation, scale = Matrix(rotation_matrix).decompose()
    return rotation


def _get_ops_from_cif(categories):
    is_petworld = False
    assembly_gen = categories["pdbx_struct_assembly_gen"]
    gen_arr = np.column_stack(
        list([assembly_gen[name].as_array() for name in assembly_gen])
    )
    dtype = [
        ("assembly_id", int),
        ("chain_id", "U10"),
        ("trans_id", int),
        ("rotation", float, 4),  # quaternion form rotations
        ("translation", float, 3),
    ]
    ops = categories["pdbx_struct_oper_list"]
    ok_names = [
        "matrix[1][1]",
        "matrix[1][2]",
        "matrix[1][3]",
        "matrix[2][1]",
        "matrix[2][2]",
        "matrix[2][3]",
        "matrix[3][1]",
        "matrix[3][2]",
        "matrix[3][3]",
        "vector[1]",
        "vector[2]",
        "vector[3]",
    ]
    # test if petworld
    if "PDB_model_num" in assembly_gen:
        print("PetWorld!")
        is_petworld = True
    op_ids = ops["id"].as_array()
    struct_ops = np.column_stack(
        list([ops[name].as_array().reshape((ops.row_count, 1))
              for name in ok_names])
    )
    rotations = np.array(
        list([rotation_from_matrix(x[0:9].reshape((3, 3)))
              for x in struct_ops])
    )
    translations = struct_ops[:, 9:12]

    gen_list = []
    for i, gen in enumerate(gen_arr):
        ids = []
        if "-" in gen[1]:
            if "," in gen[1]:
                for gexpr in gen[1].split(","):
                    if "-" in gexpr:
                        start, end = [int(x)
                                      for x in gexpr.strip("()").split("-")]
                        ids.extend((np.array(range(start, end + 1))).tolist())
                    else:
                        ids.append(int(gexpr.strip("()")))
            else:
                start, end = [int(x) for x in gen[1].strip("()").split("-")]
                ids.extend((np.array(range(start, end + 1))).tolist())
        else:
            ids = np.array([int(x)
                            for x in gen[1].strip("()").split(",")]).tolist()
        real_ids = np.nonzero(np.in1d(op_ids, [str(num) for num in ids]))[0]
        chains = np.array(gen[2].strip(" ").split(","))
        if is_petworld:
            # all chain of the model receive theses transformation
            chains = np.array([gen[3]])
        arr = np.zeros(chains.size * len(real_ids), dtype=dtype)
        arr["chain_id"] = np.tile(chains, len(real_ids))
        mask = np.repeat(np.array(real_ids), len(chains))
        if len(mask) == 0:
            print("chains are ", chains, real_ids, mask)
        try:
            arr["trans_id"] = gen[3]
        except IndexError:
            pass
        arr["rotation"] = rotations[mask, :]
        arr["translation"] = translations[mask, :]
        gen_list.append(arr)
    return np.concatenate(gen_list)