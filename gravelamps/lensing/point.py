'''
Point Mass Lensing Functions

These functions perform calculations for the isolated point mass lensing model.
Module also specifies that the C++ backend executable is pointlens

Written by Mick Wright 2022
'''

import ctypes
import os

import numpy as np

from gravelamps.core.conversion import (lens_mass_to_redshifted_lens_mass,
                                   frequency_to_dimensionless_frequency)
from gravelamps.core.gravelog import gravelogger
from .generic import get_additional_arguments

#The following loads the DLL containing the C++ functions for the direct implementations for
#geometric optics runs. It then sets the argument and result types accordingly
_cdll = ctypes.CDLL(f"{os.path.expanduser('~')}/.local/lib/libpoint.so")

_cdll.PyAmplificationFactorGeometric.argtypes = (ctypes.c_double, ctypes.c_double)
_cdll.PyAmplificationFactorGeometric.restype = ctypes.POINTER(ctypes.c_double)

_cdll.GenerateLensData.argtypes = (ctypes.c_char_p,
                                   ctypes.c_char_p,
                                   ctypes.c_char_p,
                                   ctypes.c_char_p,
                                   ctypes.c_int,
                                   ctypes.c_int)
_cdll.GenerateLensData.restype = ctypes.c_int

#Additional arguments necessary for the running of the executable in addition to the files for the
#interpolator construction
_additional_arguments = ["arithmetic_precision", "geometric_optics_frequency"]
_additional_argument_types = [int, int]

#Parameters for the model
_lens_parameters = ["lens_mass", "lens_fractional_distance", "source_position"]

def amplification_factor(dimensionless_frequency_array, source_position):
    '''
    Input:
        dimensionless_frequency_array - array of dimensionless form of the frequncies being
                                        amplified
        source_position - dimensionless displacement from the optical axis

    Output:
        amplification_array - complex values of the amplification factor for each dimensionles
                              frequency

    Function uses the C++ functions within libpoint to calculate the amplification factor in the
    geometric optics approximation for the isolated point mass model for the given dimensionless
    frequencies and source position.
    '''

    amplification_array = np.empty(len(dimensionless_frequency_array), dtype=complex)

    for idx, dimensionless_frequency in enumerate(dimensionless_frequency_array):
        c_result = _cdll.PyAmplificationFactorGeometric(ctypes.c_double(dimensionless_frequency),
                                                        ctypes.c_double(source_position))
        amplification_array[idx] = complex(c_result[0], c_result[1])
        _cdll.destroyObj(c_result)

    return amplification_array

def generate_interpolator_data(config,
                               args,
                               file_dict):
    '''
    Input:
        config - INI configuration parser
        args - Commandline argumnets passed to the program
        file_dict - dictionary of files containing dimensionless frequency and source position
                    values over which to generate the interpolator, followed by the corresponding
                    files containing the real and imaginary amplification factor values to use as
                    the interpolating data

    Function uses the C++ backend function within libpoint to generate the amplification factor
    data files from the input dimensionless frequency and source position files
    '''

    additional_arguments = get_additional_arguments(config, args,
                                                    _additional_arguments,
                                                    _additional_argument_types)

    gravelogger.info("Generating Lens Interpolator Data")
    _cdll.GenerateLensData(ctypes.c_char_p(file_dict["dimensionless_frequency"].encode("utf-8")),
                           ctypes.c_char_p(file_dict["source_position"].encode("utf-8")),
                           ctypes.c_char_p(file_dict["amplification_factor_real"].encode("utf-8")),
                           ctypes.c_char_p(file_dict["amplification_factor_imag"].encode("utf-8")),
                           ctypes.c_int(additional_arguments[0]),
                           ctypes.c_int(additional_arguments[1]))
    gravelogger.info("Lens Interpolator Data Generated")
