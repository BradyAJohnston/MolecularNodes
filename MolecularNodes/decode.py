import numpy as np
import biotite.structure as struc
import biotite.structure.io.mmtf as mmtf

types = (
    np.int8, np.int16, np.int32, 
    np.uint8, np.uint16, np.uint32, 
    np.float32, np.float64
    )
def _byte_array(bytes, enc):
    type = int(enc['type'] - 1)
    return np.frombuffer(bytes, dtype = types[type])

def _fixed_point(data, enc):
    factor = enc['factor']
    numbers = data / factor
    return numbers.astype(float)

def _interval_quant(data, enc):
    min = enc['min'] 
    max = enc['max'] 
    numsteps = enc['numsteps'] 
    srcType = enc['srcType']
    
    types = (np.float32, np.float64)
    interval = max - min / numsteps
    nums = data * interval + min
    return np.astype(nums, dtype = types[srcType])

def _run_length(data, enc):
    srcSize = enc['srcSize']
    srcType = enc['srcType']
    dtype = types[int(srcType - 1)]
    data = data.astype(dtype)
    n_entries = int(len(data) / 2)
    numbers = data.reshape((n_entries, 2))
    arr = np.repeat(numbers[:, 0], numbers[:, 1])
    return arr

def _delta(array, enc):
    origin = enc['origin']
    return np.cumsum(array) + origin

def _integer_packing(data, enc):
    # TODO currently not looking at signs - don't know how to handle this
    byteCount = int(enc['byteCount'] - 1)
    srcSize = enc['srcSize'] 
    isUnsigned = enc['isUnsigned']
    
    if isUnsigned:
        dtype = (np.uint8, np.uint16)[byteCount - 1]
    else:
        dtype = (np.int8, np.int16)[byteCount - 1]
    
    if isinstance(data, bytes):
        data = np.frombuffer(data, dtype = dtype)
    
    if len(data) > srcSize:
        if isUnsigned:
            max_sizes = (255, 65535)
        else:
            max_sizes = (127, 32767)
        max = max_sizes[byteCount]
        
        arr = np.zeros(srcSize, dtype = np.int32)
        
        idx = 0
        counter = 0
        ticker = 0
        for i, x in enumerate(data):
            if ticker < srcSize:
                if x == max:
                    counter += x
                else:
                    arr[idx] = x + counter
                    idx += 1
                    ticker += 1
                    counter = 0
                    # ticker = 0
    else:
        arr = data.astype(np.int32)
    
    return arr

def _sub_from_string(string, offsets):
    lower = 0
    substrings = []
    for upper in offsets[1:]: # drop the first range which is always zero
        substrings.append(string[lower:upper])
        lower = upper
    return np.array(substrings)

def _decode_numeric(arr, encodings):
    for enc in reversed(encodings):
        decoder = decoder_numeric[enc['kind']]
        arr = decoder(arr, enc)
    return arr

def _string_array(data, enc):
    dataEncoding = enc['dataEncoding']
    stringData = enc['stringData']
    offsetEncoding = enc['offsetEncoding']
    offsets = enc['offsets']
    
    offsets    = _decode_numeric(offsets, offsetEncoding)
    idxs       = _decode_numeric(data, dataEncoding)
    substrings = _sub_from_string(stringData, offsets)
    
    idxs[idxs >= len(substrings)] = 0 # these values should get masked further up
    
    return substrings[idxs]

decoder_numeric = {
    "ByteArray" : _byte_array,
    "FixedPoint" : _fixed_point,
    "IntervalQuantization" :  _interval_quant,
    "RunLength" :  _run_length,
    "Delta" :_delta,
    "IntegerPacking" :  _integer_packing
}

def decode(arr, encodings):
    kinds = [enc['kind'] for enc in encodings]
    if kinds[0] == 'StringArray':
        return _string_array(arr, encodings[0])
    else:
        return _decode_numeric(arr, encodings)

def decode_data(column):
    data = column['data']
    arr = data['data']
    
    mask = column.get('mask')
    
    array = decode(arr, data['encoding'])
    if mask:
        arr_mask = mask['data']
        arr_mask = _decode_numeric(arr_mask, mask['encoding'])
        return array[arr_mask]
    return array

def decode_columns(columns):
    unpacked = {}
    for col in columns:
        name = col['name']
        data = decode_data(col)
        unpacked[name] = data
    return unpacked

def get_cat_names(file):
    cats = file['dataBlocks'][0]['categories']
    names = [cat['name'] for cat in cats]
    return names

def get_cat_idx(file, cat_name):
    names = get_cat_names(file)
    return np.where(np.isin(names, cat_name))[0]

def get_atom_sites(file):
    file = mmtf.MMTFFile.read(file)
    
    idx = get_cat_idx(file, '_atom_site')[0]
    
    return decode_columns(file['dataBlocks'][0]['categories'][idx]['columns'])

def construct_atom_array(dic):
    n_atoms = len(list(dic.items())[0][1])
    names = dic.keys()
    mol = struc.AtomArray(n_atoms)
    
    coords = np.column_stack([
        dic[f'Cartn_{axis}'] for axis in "xyz"
    ])
    mol.coord = coords
    mol.chain_id = dic['label_asym_id']
    mol.atom_name = dic['label_atom_id']
    mol.element = dic['type_symbol']
    mol.res_name = dic['label_comp_id']
    mol.b_factor = dic['occupancy']
    mol.id = dic['id']
    
    return mol

def read_bcif(file):
    return construct_atom_array(get_atom_sites(file))

# fl1 = 'C:\\Users\\BradyJohnston\\git\\cellpack\\bcif\\tests\\data\\1bna.bcif'
# fl2 = "C:\\Users\\BradyJohnston\\git\\cellpack\\models\\synvesicle_2-no_bonds.bcif"


# print(construct_atom_array(get_atom_sites(fl1))[0:30])

# print(
#     construct_atom_array(get_atom_sites(fl2))[0:30]
# )
