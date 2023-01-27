"""Singular Isothermal Sphere (SIS) Lensing Functions

Following are functions performing calculations for the singular isothermal sphere lens mass
density profile model. The module backend is based in libsis.

Written by Mick Wright 2022

Globals
-------
_cdll : ctypes.CDLL
    Library of C++ functions necessary for calculations
_additional_arguments : list of str
    Additional arguments required to construct interpolator data
_additional_argument_types : list of types
    Types of the argumnets that are given above
_lens_parameters : list of str
    Parameters used for the model

Routines
--------
amplification_factor
    Calculates geometric optics amplification factor
generate_interpolator_data
    Generates the amplification factor data files for use in interpolator generation
"""

import ctypes
import os

from importlib.resources import files

import numpy as np

from gravelamps.core.gravelog import gravelogger
from .generic import get_additional_arguments

#The following loads in the DLL containing the C++ functions for the direct implementations for
#geometric optics runs. It then sets the argument and result types accordingly

_cdll = ctypes.CDLL(files('gravelamps.model.lib').joinpath('libsis.so'))

_cdll.PyAmplificationFactorGeometric.argtypes = (ctypes.c_double, ctypes.c_double)
_cdll.PyAmplificationFactorGeometric.restype = ctypes.POINTER(ctypes.c_double)

_cdll.GenerateLensData.argtypes = (ctypes.c_char_p,
                                   ctypes.c_char_p,
                                   ctypes.c_char_p,
                                   ctypes.c_char_p,
                                   ctypes.c_int,
                                   ctypes.c_int,
                                   ctypes.c_int)
_cdll.GenerateLensData.restype = ctypes.c_int

#Additional arguments necessary for the running of the executable in addition to the files for the
#interpolator construction
_additional_arguments = ["sis_summation_upper_limit",
                         "arithmetic_precision",
                         "geometric_optics_frequency"]
_additional_argument_types = [int, int, int]

#Parameters for the model
_lens_parameters = ["lens_mass", "lens_fractional_distance", "source_position"]

def amplification_factor(dimensionless_frequency_array, source_position):
    """
    Calculates geometric optics amplification factor.

    This calculation is done using C++ function `PyAmplificationFactorGeometric` within `libsis`
    for the given dimensionless frequency and source position.

    Parameters
    ----------
    dimensionless_frequency_array : Array of floats
        Dimensionless form of the frequencies of interest
    source_position : float
        Dimensionless displacement from the optical axis

    Returns
    -------
    amplification_array : Array of complex
        Amplification factor to the signal

    """

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
    """
    Generates the amplification factor data files for use in interpolator generation.

    This is done via the C++ `GenerateLensData` function within `libsis`. It will read in the
    specified grid files and fill the data files with the appropriate values of the amplification
    factor. This can be done in wave and geometric optics.

    Parameters
    ----------
    config : configparser.ConfigParser
        Object containing settings from user INI file
    args : argparse.Namespace
        Object containing commandline arguments to program
    file_dict : dict
        Contains paths to the interpolator grid and data files to fill

    """

    additional_arguments = get_additional_arguments(config, args,
                                                    _additional_arguments,
                                                    _additional_argument_types)

    gravelogger.info("Generating Lens Interpolator Data")
    _cdll.GenerateLensData(ctypes.c_char_p(file_dict["dimensionless_frequency"].encode("utf-8")),
                           ctypes.c_char_p(file_dict["source_position"].encode("utf-8")),
                           ctypes.c_char_p(file_dict["amplification_factor_real"].encode("utf-8")),
                           ctypes.c_char_p(file_dict["amplification_factor_imag"].encode("utf-8")),
                           ctypes.c_int(additional_arguments[0]),
                           ctypes.c_int(additional_arguments[1]),
                           ctypes.c_int(additional_arguments[2]))
    gravelogger.info("Lens Interpolator Data Generated")
