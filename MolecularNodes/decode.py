import numpy as np


def _byte_array(bytes, enc):
    type = int(enc['type'])
    types = (
        np.int8, np.int16, np.int32, 
        np.uint8, np.uint16, np.uint32, 
        np.float32, np.float64
        )
    return np.frombuffer(bytes, dtype = types[type])

def _fixed_point(data, enc):
    factor = enc['factor']
    srcType = enc['srcType']
    types = (np.float32, np.float64)
    numbers = data / factor
    return numbers.astype(types[srcType])

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
    n_entries = int(len(data) / 2)
    numbers = data.reshape((n_entries, 2))
    return np.repeat(numbers[0], numbers[1])

def _delta(data, enc):
    origin = enc['origin']
    return np.cumsum(array) + origin

def _integer_packing(data, enc):
    byteCount = enc['byteCount'] 
    srcSize = enc['srcSize'] 
    isUnsigned = enc['isUnsigned'] 
    dtype = (np.int8, np.int16)[byteCount]
    
    return data.astype(np.int32)
    
    
    
    
    
decoder_kind = {
    "ByteArray" : _byte_array,
    "FixedPoint" : _fixed_point,
    "IntervalQuantization" :  _interval_quant,
    "RunLength" :  _run_length,
    "Delta" :_delta,
    "IntegerPacking" :  _integer_packing
    # "StringArray" : 
}
    

dat = {
    'name': 'assembly_id',
    'data': {
        'encoding': [
            {'kind': 'RunLength', 'srcType': 3, 'srcSize': 24},
            {'kind': 'IntegerPacking','byteCount': 1,'isUnsigned': True,'srcSize': 2},
            {'kind': 'ByteArray', 'type': 4}
            ],
        'data': b'\x01\x18'
        },
    'mask': None
    }

# def run_decode(dat, enc):
#     kind = enc['kind']
#     return types[]

def decode(dic):
    data = dic['data']
    dat = data['data']
    print(dat)
    # encodings = data['encoding']
    for enc in reversed(data['encoding']):
        decoder = decoder_kind[enc['kind']]
        print(dat)
        dat = decoder(dat, enc)
        print(dat)
    return dat

print(decode(dat))