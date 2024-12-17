import numpy as np
from mathutils import Matrix
from typing import Any, Dict, List, Optional, TypedDict, Union
from biotite.structure import AtomArray


class BCIF:
    def __init__(self, file_path):
        # super().__init__()
        self.file_path = file_path
        self.file = self.read()
        self.array = _atom_array_from_bcif(self.file)
        self._transforms_data = _get_ops_from_bcif(self.file)
        self.n_models = 1
        self.n_atoms = self.array.shape
        self.chain_ids = self._chain_ids()

    def read(self):
        # if isinstance(self.file_path, BytesIO):
        #     open_bcif = self.file_path.getvalue()
        # else:
        with open(self.file_path, "rb") as data:
            open_bcif = loads(data.read())

        return open_bcif

    def assemblies(self, as_array=True):
        return self._transforms_data

    def _chain_ids(self, as_int=False) -> np.ndarray:
        if as_int:
            return np.unique(self.array.chain_id, return_inverse=True)[1]
        return np.unique(self.array.chain_id)


def _atom_array_from_bcif(open_bcif):
    categories = open_bcif.data_blocks[0]

    # check if a petworld CellPack model or not
    is_petworld = False
    if "PDB_model_num" in categories["pdbx_struct_assembly_gen"].field_names:
        print("PetWorld!")
        is_petworld = True

    atom_site = categories["atom_site"]
    n_atoms = atom_site.row_count

    # Initialise the atom array that will contain all of the data for the atoms
    # in the bcif file. TODO support multi-model bcif files
    # we first pull out the coordinates as they are from 3 different fields, but all
    # other fields should be single self-contained fields
    mol = AtomArray(n_atoms)
    coord_field_names = [f"Cartn_{axis}" for axis in "xyz"]
    mol.coord = np.hstack(
        list(
            [
                np.array(atom_site[column]).reshape((n_atoms, 1))
                for column in coord_field_names
            ]
        )
    )

    # the list of current
    atom_site_lookup = {
        # have to make sure the chain_id ends up being the same as the space operatore
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
        # annotations[0][1] = 'pdbx_PDB_model_num'
        atom_site_lookup.pop("label_asym_id")
        atom_site_lookup["pdbx_PDB_model_num"] = "chain_id"

    for name in atom_site.field_names:
        # the coordinates have already been extracted so we can skip over those field names
        if name in coord_field_names:
            continue

        # numpy does a pretty good job of guessing the data types from the fields
        data = np.array(atom_site[name])

        # if a specific name for an annotation is already specified earlier, we can
        # use that to ensure consitency. All other fields are also still added as we
        # may as well do so, in case we want any extra data
        annotation_name = atom_site_lookup.get(name)
        if not annotation_name:
            annotation_name = name

        # TODO this could be expanded to capture fields that are entirely '' and drop them
        # or fill them with 0s
        if annotation_name == "res_id" and data[0] == "":
            data = np.array([0 if x == "" else x for x in data])

        mol.set_annotation(annotation_name, data)

    return mol


def rotation_from_matrix(matrix):
    rotation_matrix = np.identity(4, dtype=float)
    rotation_matrix[:3, :3] = matrix
    translation, rotation, scale = Matrix(rotation_matrix).decompose()
    return rotation


def _get_ops_from_bcif(open_bcif):
    is_petworld = False
    cats = open_bcif.data_blocks[0]
    assembly_gen = cats["pdbx_struct_assembly_gen"]
    gen_arr = np.column_stack(
        list([assembly_gen[name] for name in assembly_gen.field_names])
    )
    dtype = [
        ("assembly_id", int),
        ("chain_id", "U10"),
        ("trans_id", int),
        ("rotation", float, 4),  # quaternion form rotations
        ("translation", float, 3),
    ]
    ops = cats["pdbx_struct_oper_list"]
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
    if "PDB_model_num" in assembly_gen.field_names:
        print("PetWorld!")
        is_petworld = True
    op_ids = np.array(ops["id"])
    struct_ops = np.column_stack(
        list([np.array(ops[name]).reshape((ops.row_count, 1)) for name in ok_names])
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
        real_ids = np.nonzero(np.in1d(op_ids, ids))[0]
        chains = np.array(gen[2].strip(" ").split(","))
        if is_petworld:
            # all chain of the model receive theses transformation
            chains = np.array([gen[3]])
        arr = np.zeros(chains.size * len(real_ids), dtype=dtype)
        arr["chain_id"] = np.tile(chains, len(real_ids))
        mask = np.repeat(np.array(real_ids), len(chains))
        try:
            arr["trans_id"] = gen[3]
        except IndexError as e:
            print(e)
        arr["rotation"] = rotations[mask, :]
        arr["translation"] = translations[mask, :]
        gen_list.append(arr)
    return np.concatenate(gen_list)


# This BinaryCIF implementation was taken from here: https://gist.github.com/dsehnal/b06f5555fa9145da69fe69abfeab6eaf

# BinaryCIF Parser
# Copyright (c) 2021 David Sehnal <david.sehnal@gmail.com>, licensed under MIT.
#
# Resources:
# - https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008247
# - https://github.com/molstar/BinaryCIF & https://github.com/molstar/BinaryCIF/blob/master/encoding.md
#
# Implementation based on Mol*:
# - https://github.com/molstar/molstar/blob/master/src/mol-io/common/binary-cif/encoding.ts
# - https://github.com/molstar/molstar/blob/master/src/mol-io/common/binary-cif/decoder.ts
# - https://github.com/molstar/molstar/blob/master/src/mol-io/reader/cif/binary/parser.ts


class EncodingBase(TypedDict):
    kind: str


class EncodedData(TypedDict):
    encoding: List[EncodingBase]
    data: bytes


class EncodedColumn(TypedDict):
    name: str
    data: EncodedData
    mask: Optional[EncodedData]


class EncodedCategory(TypedDict):
    name: str
    rowCount: int
    columns: List[EncodedColumn]


class EncodedDataBlock(TypedDict):
    header: str
    categories: List[EncodedCategory]


class EncodedFile(TypedDict):
    version: str
    encoder: str
    dataBlocks: List[EncodedDataBlock]


def _decode(encoded_data: EncodedData) -> Union[np.ndarray, List[str]]:
    result = encoded_data["data"]
    for encoding in encoded_data["encoding"][::-1]:
        if encoding["kind"] in _decoders:
            result = _decoders[encoding["kind"]](result, encoding)  # type: ignore
        else:
            raise ValueError(f"Unsupported encoding '{encoding['kind']}'")

    return result  # type: ignore


class DataTypes:
    Int8 = 1
    Int16 = 2
    Int32 = 3
    Uint8 = 4
    Uint16 = 5
    Uint32 = 6
    Float32 = 32
    Float64 = 33


_dtypes = {
    DataTypes.Int8: "i1",
    DataTypes.Int16: "i2",
    DataTypes.Int32: "i4",
    DataTypes.Uint8: "u1",
    DataTypes.Uint16: "u2",
    DataTypes.Uint32: "u4",
    DataTypes.Float32: "f4",
    DataTypes.Float64: "f8",
}


def _get_dtype(type: int) -> str:
    if type in _dtypes:
        return _dtypes[type]

    raise ValueError(f"Unsupported data type '{type}'")


class ByteArrayEncoding(EncodingBase):
    type: int


class FixedPointEncoding(EncodingBase):
    factor: float
    srcType: int


class IntervalQuantizationEncoding(EncodingBase):
    min: float
    max: float
    numSteps: int
    srcType: int


class RunLengthEncoding(EncodingBase):
    srcType: int
    srcSize: int


class DeltaEncoding(EncodingBase):
    origin: int
    srcType: int


class IntegerPackingEncoding(EncodingBase):
    byteCount: int
    isUnsigned: bool
    srcSize: int


class StringArrayEncoding(EncodingBase):
    dataEncoding: List[EncodingBase]
    stringData: str
    offsetEncoding: List[EncodingBase]
    offsets: bytes


def _decode_byte_array(data: bytes, encoding: ByteArrayEncoding) -> np.ndarray:
    return np.frombuffer(data, dtype="<" + _get_dtype(encoding["type"]))


def _decode_fixed_point(data: np.ndarray, encoding: FixedPointEncoding) -> np.ndarray:
    return np.array(data, dtype=_get_dtype(encoding["srcType"])) / encoding["factor"]


def _decode_interval_quantization(
    data: np.ndarray, encoding: IntervalQuantizationEncoding
) -> np.ndarray:
    delta = (encoding["max"] - encoding["min"]) / (encoding["numSteps"] - 1)
    return (
        np.array(data, dtype=_get_dtype(encoding["srcType"])) * delta + encoding["min"]
    )


def _decode_run_length(data: np.ndarray, encoding: RunLengthEncoding) -> np.ndarray:
    return np.repeat(
        np.array(data[::2], dtype=_get_dtype(encoding["srcType"])), repeats=data[1::2]
    )


def _decode_delta(data: np.ndarray, encoding: DeltaEncoding) -> np.ndarray:
    result = np.array(data, dtype=_get_dtype(encoding["srcType"]))
    if encoding["origin"]:
        result[0] += encoding["origin"]
    return np.cumsum(result, out=result)


def _decode_integer_packing_signed(
    data: np.ndarray, encoding: IntegerPackingEncoding
) -> np.ndarray:
    upper_limit = 0x7F if encoding["byteCount"] == 1 else 0x7FFF
    lower_limit = -upper_limit - 1
    n = len(data)
    output = np.zeros(encoding["srcSize"], dtype="i4")
    i = 0
    j = 0
    while i < n:
        value = 0
        t = data[i]
        while t == upper_limit or t == lower_limit:
            value += t
            i += 1
            t = data[i]
        value += t
        output[j] = value
        i += 1
        j += 1
    return output


def _decode_integer_packing_unsigned(
    data: np.ndarray, encoding: IntegerPackingEncoding
) -> np.ndarray:
    upper_limit = 0xFF if encoding["byteCount"] == 1 else 0xFFFF
    n = len(data)
    output = np.zeros(encoding["srcSize"], dtype="i4")
    i = 0
    j = 0
    while i < n:
        value = 0
        t = data[i]
        while t == upper_limit:
            value += t
            i += 1
            t = data[i]
        value += t
        output[j] = value
        i += 1
        j += 1
    return output


def _decode_integer_packing(
    data: np.ndarray, encoding: IntegerPackingEncoding
) -> np.ndarray:
    if len(data) == encoding["srcSize"]:
        return data
    if encoding["isUnsigned"]:
        return _decode_integer_packing_unsigned(data, encoding)
    else:
        return _decode_integer_packing_signed(data, encoding)


def _decode_string_array(data: bytes, encoding: StringArrayEncoding) -> List[str]:
    offsets = _decode(
        EncodedData(encoding=encoding["offsetEncoding"], data=encoding["offsets"])
    )
    indices = _decode(EncodedData(encoding=encoding["dataEncoding"], data=data))

    str = encoding["stringData"]
    strings = [""]
    for i in range(1, len(offsets)):
        strings.append(str[offsets[i - 1] : offsets[i]])  # type: ignore

    return [strings[i + 1] for i in indices]  # type: ignore


_decoders = {
    "ByteArray": _decode_byte_array,
    "FixedPoint": _decode_fixed_point,
    "IntervalQuantization": _decode_interval_quantization,
    "RunLength": _decode_run_length,
    "Delta": _decode_delta,
    "IntegerPacking": _decode_integer_packing,
    "StringArray": _decode_string_array,
}


###############################################################################


class CifValueKind:
    Present = 0
    # Expressed in CIF as `.`
    NotPresent = 1
    # Expressed in CIF as `?`
    Unknown = 2


class CifField:
    def __getitem__(self, idx: int) -> Union[str, float, int, None]:
        # if self._value_kinds and self._value_kinds[idx]:
        # return None
        return self._values[idx]

    def __len__(self):
        return self.row_count

    @property
    def values(self):
        """
        A numpy array of numbers or a list of strings.
        """
        return self._values

    @property
    def value_kinds(self):
        """
        value_kinds represent the presence or absence of particular "CIF value".
        - If the mask is not set, every value is present:
            - 0 = Value is present
            - 1 = . = value not specified
            - 2 = ? = value unknown
        """
        return self._value_kinds

    def __init__(
        self,
        name: str,
        values: Union[np.ndarray, List[str]],
        value_kinds: Optional[np.ndarray],
    ):
        self.name = name
        self._values = values
        self._value_kinds = value_kinds
        self.row_count = len(values)


class CifCategory:
    def __getattr__(self, name: str) -> Any:
        return self[name]

    def __getitem__(self, name: str) -> Optional[CifField]:
        if name not in self._field_cache:
            return None

        if not self._field_cache[name]:
            self._field_cache[name] = _decode_column(self._columns[name])

        return self._field_cache[name]

    def __contains__(self, key: str):
        return key in self._columns

    def __init__(self, category: EncodedCategory, lazy: bool):
        self.field_names = [c["name"] for c in category["columns"]]
        self._field_cache = {
            c["name"]: None if lazy else _decode_column(c) for c in category["columns"]
        }
        self._columns: Dict[str, EncodedColumn] = {
            c["name"]: c for c in category["columns"]
        }
        self.row_count = category["rowCount"]
        self.name = category["name"][1:]


class CifDataBlock:
    def __getattr__(self, name: str) -> Any:
        return self.categories[name]

    def __getitem__(self, name: str) -> CifCategory:
        return self.categories[name]

    def __contains__(self, key: str):
        return key in self.categories

    def __init__(self, header: str, categories: Dict[str, CifCategory]):
        self.header = header
        self.categories = categories


class CifFile:
    def __getitem__(self, index_or_name: Union[int, str]):
        """
        Access a data block by index or header (case sensitive)
        """
        if isinstance(index_or_name, str):
            return (
                self._block_map[index_or_name]
                if index_or_name in self._block_map
                else None
            )
        else:
            return (
                self.data_blocks[index_or_name]
                if index_or_name < len(self.data_blocks)
                else None
            )

    def __len__(self):
        return len(self.data_blocks)

    def __contains__(self, key: str):
        return key in self._block_map

    def __init__(self, data_blocks: List[CifDataBlock]):
        self.data_blocks = data_blocks
        self._block_map = {b.header: b for b in data_blocks}


def _decode_column(column: EncodedColumn) -> CifField:
    values = _decode(column["data"])
    value_kinds = _decode(column["mask"]) if column["mask"] else None  # type: ignore
    # type: ignore
    return CifField(name=column["name"], values=values, value_kinds=value_kinds)


def loads(data: Union[bytes, EncodedFile], lazy=True) -> CifFile:
    """
    - data: msgpack encoded blob or EncodedFile object
    - lazy:
        - True: individual columns are decoded only when accessed
        - False: decode all columns immediately
    """
    import msgpack

    file: EncodedFile = (
        data if isinstance(data, dict) and "dataBlocks" in data else msgpack.loads(data)
    )  # type: ignore

    data_blocks = [
        CifDataBlock(
            header=block["header"],
            categories={
                cat["name"][1:]: CifCategory(category=cat, lazy=lazy)
                for cat in block["categories"]
            },
        )
        for block in file["dataBlocks"]
    ]

    return CifFile(data_blocks=data_blocks)
