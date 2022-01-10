import numpy as np
import os
import ctypes

_sis = ctypes.CDLL(os.path.abspath("libsis.so"))
_sis.AFGRealOnly.argtypes = (ctypes.c_double, ctypes.c_double) 
_sis.AFGRealOnly.restype = ctypes.POINTER(ctypes.c_double) 

def amplification_factor_geometric(dimensionless_frequency, source_position):
    global _sis

    result = _sis.AFGRealOnly(ctypes.c_double(dimensionless_frequency), ctypes.c_double(source_position))
    py_result = complex(result[0], result[1])

    _sis.destroyObj(result) 

    return result 

amplification_factor_geometric(1,1)
