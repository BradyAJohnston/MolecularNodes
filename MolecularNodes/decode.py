import numpy as np


def _byte_array(bytes, enc):
    type = int(enc['type'] - 1)
    types = (
        np.int8, np.int16, np.int32, 
        np.uint8, np.uint16, np.uint32, 
        np.float32, np.float64
        )
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
    is_odd = len(data) % 2 != 0
    if is_odd:
        data = data[:-1]
    n_entries = int(len(data) / 2)
    numbers = data.reshape((n_entries, 2))
    return np.repeat(numbers[:, 0], numbers[:, 1])

def _delta(array, enc):
    origin = enc['origin']
    return np.cumsum(array) + origin

def _integer_packing(data, enc):
    # TODO currently not looking at signs - don't know how to handle this
    byteCount = int(enc['byteCount'] - 1)
    srcSize = enc['srcSize'] 
    isUnsigned = enc['isUnsigned'] 
    dtype = (np.int8, np.int16)[byteCount]
    if isinstance(data, bytes):
        data = np.frombuffer(data, dtype = dtype)
    
    if len(data) > srcSize:
        if byteCount == 0:
            if isUnsigned:
                max = 255
            else:
                max = 127
        elif byteCount == 1:
            if isUnsigned:
                max = 65535
            else:
                max = 32768
        
        arr = np.zeros(srcSize, dtype = np.int32)
        
        idx = 0
        counter = 0
        for i, x in enumerate(data):
            if x == 255:
                counter += x
                next
            else:
                arr[idx] = x + counter
                idx += 1
                counter = 0
    else:
        arr = data.astype(np.int32)
    
    return arr

def _sub_from_string(string, offsets):
    lower = 0
    substrings = []
    for upper in offsets:
        substrings.append(string[lower:upper])
        lower = upper
    return np.array(substrings)

def _decode_numeric(arr, encodings):
    for enc in reversed(encodings):
        decoder = decoder_numeric[enc['kind']]
        arr = decoder(arr, enc)
        print(arr)
    return arr

def _string_array(data, enc):
    dataEncoding = enc['dataEncoding']
    stringData = enc['stringData']
    offsetEncoding = enc['offsetEncoding']
    offsets = enc['offsets']
    
    offsets = _decode_numeric(offsets, offsetEncoding)
    print(offsets)
    print(dataEncoding)
    print(data)
    idxs   = _decode_numeric(data, dataEncoding)
    print(idxs)
    substrings = _sub_from_string(stringData, offsets)
    print(substrings)
    return substrings[idxs]



decoder_numeric = {
    "ByteArray" : _byte_array,
    "FixedPoint" : _fixed_point,
    "IntervalQuantization" :  _interval_quant,
    "RunLength" :  _run_length,
    "Delta" :_delta,
    "IntegerPacking" :  _integer_packing
}

dat = {
    'name': 'auth_atom_id',
    'data': {
        'encoding': [
            {
                'kind': 'StringArray',
                'dataEncoding': [
                    {'kind': 'Delta', 'origin': 0, 'srcType': 3},
                    {'kind': 'RunLength', 'srcType': 3, 'srcSize': 566},
                    {'kind': 'IntegerPacking', 'byteCount': 1, 'isUnsigned': False, 'srcSize': 412},
                    {'kind': 'ByteArray', 'type': 1}
                ],
                'stringData': "O5'C5'C4'O4'C3'O3'C2'C1'N1C2O2N3C4N4C5C6POP1OP2N9C8N7O6N2N6O4C7O",
                'offsetEncoding': [
                    {'kind': 'Delta', 'origin': 0, 'srcType': 3},
                    {'kind': 'RunLength', 'srcType': 3, 'srcSize': 29},
                    {'kind': 'IntegerPacking', 'byteCount': 1, 'isUnsigned': True, 'srcSize': 14},
                    {'kind': 'ByteArray', 'type': 4}
                ],
                'offsets': b'\x00\x01\x03\x08\x02\x08\x01\x01\x03\x02\x02\x08\x01\x01'
            }
        ],
        'data': b'\x00\x01\x01\x12\xee\x01\x01\x07\x0c\x01\x01\x02\xf9\x01\x01\x01\x07\x01\xf2\x01\x01\x01\x0e\x01\xf4\x01\x01\x01\x04\x01\x01\x02\xee\x01\x01\x12\xee\x01\x01\x07\x0c\x01\x01\x02\xf9\x01\x01\x01\x07\x01\xf2\x01\x01\x01\x0e\x01\xf4\x01\x01\x01\x04\x01\x01\x02\xee\x01\x01\x07\x0c\x01\x01\x02\xf9\x01\x01\x01\t\x01\xf0\x01\x01\x01\x02\x01\x01\x01\x04\x01\x01\x02\xee\x01\x01\x07\x0c\x01\x01\x02\xf9\x01\x01\x01\t\x01\xf0\x01\x01\x01\x02\x01\x01\x01\x04\x01\x01\x02\xee\x01\x01\x0c\r\x01\xf5\x01\x0c\x01\xf5\x01\x01\x03\xee\x01\x01\x0c\r\x01\xf5\x01\x0c\x01\xf5\x01\x01\x03\xee\x01\x01\x12\xee\x01\x01\x07\x0c\x01\x01\x02\xf9\x01\x01\x01\x07\x01\xf2\x01\x01\x01\x0e\x01\xf4\x01\x01\x01\x04\x01\x01\x02\xee\x01\x01\x12\xee\x01\x01\x07\x0c\x01\x01\x02\xf9\x01\x01\x01\x07\x01\xf2\x01\x01\x01\x0e\x01\xf4\x01\x01\x01\xf4\x01\x01\x12\xee\x01\x01\x07\x0c\x01\x01\x02\xf9\x01\x01\x01\x07\x01\xf2\x01\x01\x01\x0e\x01\xf4\x01\x01\x01\x04\x01\x01\x02\xee\x01\x01\x12\xee\x01\x01\x07\x0c\x01\x01\x02\xf9\x01\x01\x01\x07\x01\xf2\x01\x01\x01\x0e\x01\xf4\x01\x01\x01\x04\x01\x01\x02\xee\x01\x01\x07\x0c\x01\x01\x02\xf9\x01\x01\x01\t\x01\xf0\x01\x01\x01\x02\x01\x01\x01\x04\x01\x01\x02\xee\x01\x01\x07\x0c\x01\x01\x02\xf9\x01\x01\x01\t\x01\xf0\x01\x01\x01\x02\x01\x01\x01\x04\x01\x01\x02\xee\x01\x01\x0c\r\x01\xf5\x01\x0c\x01\xf5\x01\x01\x03\xee\x01\x01\x0c\r\x01\xf5\x01\x0c\x01\xf5\x01\x01\x03\xee\x01\x01\x12\xee\x01\x01\x07\x0c\x01\x01\x02\xf9\x01\x01\x01\x07\x01\xf2\x01\x01\x01\x0e\x01\xf4\x01\x01\x01\x04\x01\x01\x02\xee\x01\x01\x12\xee\x01\x01\x07\x0c\x01\x01\x02\xf9\x01\x01\x01\x07\x01\xf2\x01\x01\x01\x0e\x01\xf4\x01\x01\x01\x0f\x01\x00O'
    },
    'mask': None
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



dat = {
    'encoding': [
        {'kind': 'RunLength', 'srcType': 3, 'srcSize': 566}, 
        {'kind': 'IntegerPacking', 'byteCount': 1, 'isUnsigned': True, 'srcSize': 4}, 
        {'kind': 'ByteArray', 'type': 4}
        ], 
    'data': b'\x00\xff\xe7\x01P'
    }

print(
    _decode_numeric(dat['data'], dat['encoding'])
)



def decode_columns(columns):
    unpacked = {}
    for col in columns:
        name = col['name']
        data = decode_data(col)
        unpacked[name] = data
        print(unpacked)
    return unpacked

def get_cat_names(file):
    cats = file['dataBlocks'][0]['categories']
    names = [cat['name'] for cat in cats]
    return names

def get_cat_idx(file, cat_name):
    names = get_cat_names(file)
    return np.where(np.isin(names, cat_name))[0]

def get_atom_sites(file):
    import biotite.structure.io.mmtf as mmtf
    file = mmtf.MMTFFile.read(file)
    
    idx = get_cat_idx(file, '_atom_site')[0]
    print(idx)
    
    return decode_columns(file['dataBlocks'][0]['categories'][idx]['columns'])

print(
    get_atom_sites(
        'C:\\Users\\BradyJohnston\\git\\cellpack\\bcif\\tests\\data\\1bna.bcif'
        )
    )
