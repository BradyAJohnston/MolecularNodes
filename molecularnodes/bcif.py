import numpy as np
import warnings


def rotation_from_matrix(matrix):
    from scipy.spatial.transform import Rotation
    with warnings.catch_warnings():
        rotation = Rotation.from_matrix(matrix).as_euler('xyz')
    return rotation


def parse(file):
    with open(file, "rb") as data:
        open_bcif = loads(data.read())
    
    mol = atom_array_from_bcif(open_bcif)
    syms = get_ops_from_bcif(open_bcif)

    return mol, syms


def get_ops_from_bcif(open_bcif):
    is_petworld = False
    cats = open_bcif.data_blocks[0]
    assembly_gen = cats['pdbx_struct_assembly_gen']
    gen_arr = np.column_stack(list([assembly_gen[name] for name in assembly_gen.field_names]))
    dtype = [
        ('assembly_id', int),
        ('chain_id',    'U10'),
        ('trans_id', int),
        ('rotation',    float, 3),
        ('translation', float, 3)
    ]
    ops = cats['pdbx_struct_oper_list']
    ok_names = [
        'matrix[1][1]',
        'matrix[1][2]',
        'matrix[1][3]',
        'matrix[2][1]',
        'matrix[2][2]',
        'matrix[2][3]',
        'matrix[3][1]',
        'matrix[3][2]',
        'matrix[3][3]',
        'vector[1]',
        'vector[2]',
        'vector[3]'
    ]
    # test if petworld
    if 'PDB_model_num' in assembly_gen.field_names:
        print('PetWorld!')
        is_petworld = True
    ninstance = len(ops['vector[1]'])
    op_ids = np.array(ops['id'])
    struct_ops = np.column_stack(list([
        np.array(ops[name]).reshape((ops.row_count, 1)) for name in ok_names
    ]))
    rotations = np.array(list([
        rotation_from_matrix(x[0:9].reshape((3, 3))) for x in struct_ops
    ]))
    translations = struct_ops[:, 9:12]

    gen_list = []
    for i, gen in enumerate(gen_arr):
        ids = []
        if "-" in gen[1]:
            if "," in gen[1]:
                for gexpr in gen[1].split(","):
                    if "-" in gexpr:
                        start, end = [int(x) for x in gexpr.strip('()').split('-')]
                        ids.extend( (np.array(range(start, end + 1))).tolist())
                    else:
                        ids.append(int(gexpr.strip('()')))
            else:
                start, end = [int(x) for x in gen[1].strip('()').split('-')]
                ids.extend( (np.array(range(start, end + 1))).tolist())
        else:
            ids = np.array([int(x) for x in gen[1].strip("()").split(",")]).tolist()
        real_ids = np.nonzero(np.in1d(op_ids,ids))[0]
        chains = np.array(gen[2].strip(' ').split(','))
        if is_petworld:
            # all chain of the model receive theses transformation
            chains = np.array([gen[3]])
        arr = np.zeros(chains.size * len(real_ids), dtype = dtype)
        arr['chain_id']    = np.tile(chains, len(real_ids))
        mask = np.repeat(np.array(real_ids), len(chains))
        try:
            arr['trans_id']    = gen[3]
        except IndexError:
            pass
        arr['rotation']    = rotations[mask, :]
        arr['translation'] = translations[mask, :]
        gen_list.append(arr)
    return np.concatenate(gen_list)


def atom_array_from_bcif(open_bcif):
    from biotite.structure import AtomArray
    is_petworld = False
    cats = open_bcif.data_blocks[0]
    assembly_gen = cats['pdbx_struct_assembly_gen']
    if 'PDB_model_num' in assembly_gen.field_names:
        print('PetWorld!')
        is_petworld = True    
    atom_site = open_bcif.data_blocks[0].categories['atom_site']
    n_atoms = atom_site.row_count
    mol = AtomArray(n_atoms)

    coords = np.hstack(list([
        np.array(atom_site[f'Cartn_{axis}']).reshape(n_atoms, 1) for axis in 'xyz'
    ]))
    mol.coord = coords

    annotations = [
        # have to be the same
        # chainid as in space operator
        ['chain_id',  'label_asym_id'],
        ['atom_name', 'label_atom_id'],
        ['res_name',  'label_comp_id'],
        ['element',   'type_symbol'],
        ['res_id',    'label_seq_id'],  # or auth
        ['b_factor',  'B_iso_or_equiv'],
        ['entity_id', 'label_entity_id'],
        ['model_id', 'pdbx_PDB_model_num']
    ]
    if is_petworld:
        annotations[0][1] = 'pdbx_PDB_model_num'
    for ann in annotations:
        dat = atom_site[ann[1]]
        print(ann)
        if dat:
            if dat[0] == '' and ann[0] == 'res_id':
                dat = np.array([dat[x] if dat[0] != '' else '0'  for x in range(len(dat))])
            mol.set_annotation(ann[0], dat)
    return mol


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

from typing import Any, Dict, List, Optional, TypedDict, Union

import msgpack
import numpy as np


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
    return CifField(name=column["name"], values=values, value_kinds=value_kinds)  # type: ignore


def loads(data: Union[bytes, EncodedFile], lazy=True) -> CifFile:
    """
    - data: msgpack encoded blob or EncodedFile object
    - lazy:
        - True: individual columns are decoded only when accessed
        - False: decode all columns immediately
    """

    file: EncodedFile = data if isinstance(data, dict) and "dataBlocks" in data else msgpack.loads(data)  # type: ignore

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
