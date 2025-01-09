from io import BytesIO
from .bcif import BCIF
from ..molecule.pdbx import PDBX, _extract_matrices
from biotite.structure.io import pdbx
from biotite import structure as struc
from biotite import InvalidFileError
from pathlib import Path
import numpy as np


# we have to override the reading of lines to strip out the preceeding whitespace
# for the petworld files
class PetworldCIFFileReader(pdbx.CIFFile):
    @classmethod
    def read(cls, file):
        """
        Read a CIF file.

        Parameters
        ----------
        file : file-like object or str
            The file to be read.
            Alternatively a file path can be supplied.

        Returns
        -------
        file_object : CIFFile
            The parsed file.
        """
        # File name
        if pdbx.cif.is_open_compatible(file):
            with open(file, "r") as f:
                text = f.read()
        # File object
        else:
            if not pdbx.cif.is_text(file):
                raise TypeError("A file opened in 'text' mode is required")
            text = file.read()
        return cls.deserialize(text)

    @staticmethod
    def deserialize(text):
        lines = [line.strip() for line in text.splitlines()]
        block_starts = []
        block_names = []
        for i, line in enumerate(lines):
            if not pdbx.cif._is_empty(line):
                data_block_name = pdbx.cif._parse_data_block_name(line)
                if data_block_name is not None:
                    block_starts.append(i)
                    block_names.append(data_block_name)
        return pdbx.CIFFile(
            pdbx.cif._create_element_dict(lines, block_names, block_starts)
        )


class CellPackReader(PDBX):
    def __init__(self, file_path):
        super().__init__(file_path)
        self._extra_annotations["asym_id"] = self._get_asym_id
        self._extra_annotations["pdb_model_num"] = self._get_pdb_model_num
        self.file_path = file_path
        self.file: pdbx.BinaryCIFFile | pdbx.CIFFile = self._read()
        self.n_molecules: int = pdbx.get_model_count(self.file)
        self.molecules: dict[str, stuc.AtomArray] = {}

    @property
    def mol_ids(self) -> np.ndarray:
        return np.unique(list(self.molecules.keys()))

    def _parse_filepath(self, file_path: Path | str | BytesIO) -> None:
        pass

    def _read(self):
        suffix = Path(self.file_path).suffix
        if suffix == ".bcif":
            return pdbx.BinaryCIFFile.read(self.file_path)
        elif suffix == ".cif":
            return PetworldCIFFileReader.read(self.file_path)
        else:
            raise ValueError(f"Invalid file format: '{suffix}")

    def _get_asym_id(self, array, file) -> np.ndarray:
        return array.chain_id

    def _get_pdb_model_num(self, array, file) -> np.ndarray:
        return (
            list(self.file.values())[0]["atom_site"]["pdbx_PDB_model_num"]
            .as_array()
            .astype(int)
        )

    def get_structure(
        self,
        extra_fields=["b_factor", "occupancy", "atom_id"],
        bonds=True,
        model: int | None = None,
    ) -> struc.AtomArray:
        array: struc.AtomArray = pdbx.get_structure(
            self.file, model=model, extra_fields=extra_fields
        )  # type: ignore
        array = self.set_extra_annotations(array, self.file)

        if not array.bonds and bonds:
            array.bonds = struc.bonds.connect_via_residue_names(
                array, inter_residue=True
            )
        return array

    @property
    def blocks(self):
        return list(self.file.values())[0]

    def get_molecules(
        self, extra_fields=["b_factor", "occupancy", "atom_id"], bonds=True
    ):
        self._is_petworld = False

        if "PDB_model_num" in self.blocks["pdbx_struct_assembly_gen"]:
            self._is_petworld = True

        try:
            array = self.get_structure(extra_fields, bonds)

            if self._is_petworld:
                array.set_annotation(
                    "chain_id", np.char.rjust(array.pdb_model_num, 4, "0")
                )

            self.molecules = {
                c: array[array.chain_id == c] for c in np.unique(array.chain_id)
            }
        except InvalidFileError as e:
            self._is_petworld = True
            for i in range(self.n_molecules):
                array = self.get_structure(extra_fields, bonds, model=i + 1)
                array.set_annotation("pdbx_PDB_model_num", np.repeat(i + 1, len(array)))
                array.chain_id = array.pdbx_PDB_model_num
                chain_name = "{}_{}".format(
                    str(i).rjust(4, "0"), str(array.chain_id[0])
                )
                self.molecules[chain_name] = array

    def _get_petworld_assemblies(self) -> np.ndarray:
        """
        Returns the assembly ids and chain ids of the PetWorld files.
        """
        assmbly_details = self.blocks["pdbx_struct_assembly_gen"]
        matrices = _extract_matrices(self.blocks["pdbx_struct_oper_list"])
        assembly_ids = assmbly_details["assembly_id"].as_array(str)
        oper_expression = assmbly_details["oper_expression"].as_array(str)
        asym_id_list = assmbly_details["asym_id_list"].as_array(str)
        PDB_model_num = assmbly_details["PDB_model_num"].as_array(int)


def _get_ops_from_cif(categories, lookup=None):
    is_petworld = False
    assembly_gen = categories["pdbx_struct_assembly_gen"]
    gen_arr = np.column_stack(
        list([assembly_gen[name].as_array() for name in assembly_gen])
    )
    dtype = [
        ("assembly_id", int),
        ("chain_id", "U10"),
        ("transform_id", int),
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
        is_petworld = True
    # operator ID can be a string
    op_ids = ops["id"].as_array(str)
    struct_ops = np.column_stack(
        list([ops[name].as_array().reshape((ops.row_count, 1)) for name in ok_names])
    )
    rotations = np.array(
        list([rotation_from_matrix(x[0:9].reshape((3, 3))) for x in struct_ops])
    )
    translations = struct_ops[:, 9:12]

    gen_list = []
    for i, gen in enumerate(gen_arr):
        ids = []
        if "-" in gen[1]:
            if "," in gen[1]:
                for gexpr in gen[1].split(","):
                    if "-" in gexpr:
                        start, end = [int(x) for x in gexpr.strip("()").split("-")]
                        ids.extend((np.array(range(start, end + 1))).tolist())
                    else:
                        ids.append(int(gexpr.strip("()")))
            else:
                start, end = [int(x) for x in gen[1].strip("()").split("-")]
                ids.extend((np.array(range(start, end + 1))).tolist())
        else:
            ids = np.array([int(x) for x in gen[1].strip("()").split(",")]).tolist()
        real_ids = np.nonzero(np.in1d(op_ids, [str(num) for num in ids]))[0]
        chains = np.array(gen[2].strip(" ").split(","))
        if is_petworld:
            # all chain of the model receive theses transformation
            chains = np.array([gen[3]])
        arr = np.zeros(chains.size * len(real_ids), dtype=dtype)
        arr["chain_id"] = np.tile(chains, len(real_ids))
        mask = np.repeat(np.array(real_ids), len(chains))
        if len(mask) == 0:
            print("no mask chains are ", chains, real_ids, mask)
        try:
            arr["assembly_id"] = gen[0]
        except IndexError:
            pass
        if is_petworld:
            arr["transform_id"] = gen[3]
        else:
            if lookup:
                arr["transform_id"] = np.array(
                    [lookup[chain] for chain in arr["chain_id"]]
                )
            else:
                arr["transform_id"] = mask
        arr["rotation"] = rotations[mask, :]
        arr["translation"] = translations[mask, :]
        gen_list.append(arr)
    return np.concatenate(gen_list)
