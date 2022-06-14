'''
Point Mass Lensing Functions

These functions perform calculations for the isolated point mass lensing model.
Module also specifies that the C++ backend executable is pointlens

Written by Mick Wright 2022
'''

import ctypes
import os

import numpy as np

#The following loads the DLL containing the C++ functions for the direct implementations for
#geometric optics runs. It then sets the argument and result types accordingly

_cdll = ctypes.CDLL(f"{os.path.expanduser('~')}/.local/lib/libpoint.so")
_cdll.AFGRealOnly.argtypes = (ctypes.c_double, ctypes.c_double)
_cdll.AFGRealOnly.restype = ctypes.POINTER(ctypes.c_double)

#Location of the executable. This executable is made by Gravelamps upon install and placed
#into the bin directory in the user's home folder
_executable = f"{os.path.expanduser('~')}/bin/pointlens"

#Additional arguments necessary for the running of the executable in addition to the files for the
#interpolator construction
_additional_arguments = ["arithmetic_precision", "geometric_optics_frequency"]

def amplification_factor(dimensionless_frequency_array, source_position):
    '''
    Input:
        dimensionless_frequency_array - values of dimensionless frequency over which to generate
                                        the amplification factor
        source_position - dimensionless displacement from the optical axis

    Output:
        amplification_array - complex values of the amplification factor for each dimensionles
                              frequency

    Function uses the C++ functions within libpoint to calculate the amplification factor in the
    geometric optics approximation for the isolated point mass model for the given dimensionless
    frequencies and source position.
    '''

    amplification_array = np.zeros(len(dimensionless_frequency_array), dtype=complex)

    for idx, dimensionless_frequency in enumerate(dimensionless_frequency_array):
        c_result = _cdll.AFGRealOnly(ctypes.c_double(dimensionless_frequency),
                                     ctypes.c_double(source_position))
        amplification_array[idx] = complex(c_result[0], c_result[1])
        _cdll.destroyObj(c_result)

    return amplification_array
