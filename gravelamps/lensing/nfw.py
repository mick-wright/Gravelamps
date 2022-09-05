'''
Navarro Frenk White (NFW) Lensing Functions

These functions perform calculations for the Navarro, Frenk, White lensing model.
Module also specifies that the C++ backend executable is nfwlens.

Written by Mick Wright 2022
'''

from collections.abc import Callable

import ctypes
import os

import numpy as np
from numpy.typing import ArrayLike
from scipy.interpolate import interp1d

from gravelamps.core.gravelog import gravelogger
from gravelamps.lensing.generic import get_additional_arguments

#The following loads in the DLL containing the C++ functions for the direct implementations for
#goemetric optics runs. It then sets the argument and result types accordingly

_cdll = ctypes.CDLL(f"{os.path.expanduser('~')}/.local/lib/libnfw.so")

_cdll.PyPhase.argtypes = (ctypes.c_double, ctypes.c_double)
_cdll.PyPhase.restype = ctypes.c_double

_cdll.PyImagePositions.argtypes = (ctypes.c_double, ctypes.c_double)
_cdll.PyImagePositions.restype = ctypes.POINTER(ctypes.c_double)

_cdll.PyAmplificationFactorGeometric.argtypes = (ctypes.c_double,
                                                 ctypes.c_double,
                                                 ctypes.c_double,
                                                 ctypes.POINTER(ctypes.c_double),
                                                 ctypes.c_double,
                                                 ctypes.c_int)
_cdll.PyAmplificationFactorGeometric.restype = ctypes.POINTER(ctypes.c_double)

_cdll.GenerateLensData.argtypes = (ctypes.c_char_p,
                                   ctypes.c_char_p,
                                   ctypes.c_char_p,
                                   ctypes.c_char_p,
                                   ctypes.c_double,
                                   ctypes.c_double,
                                   ctypes.c_int,
                                   ctypes.c_int)
_cdll.GenerateLensData.restype = ctypes.c_int

#Additional arguments necessary for the running of the executable in addition to the files for the
#interpolator construction
_additional_arguments = ["nfw_scaling_constant",
                         "nfw_integration_upper_limit",
                         "arithmetic_precision",
                         "geometric_optics_frequency"]
_additional_argument_types = [float, float, int, int]

#Parameters for the model
_lens_parameters = ["lens_mass", "lens_fractional_distance", "source_position"]

#The following are global properties of the module. These are the scaling constant of the profile
#and the critical value of source position --- that corresponding to the caustic. The scaling may
#be adjusted by the set_scaling function which will also correspondingly calculate the appropriate
#critical value
_SCALING : float
_CRITICAL_VALUE : float

#The following are interpolators that are dynamically created whenever the scaling is adjusted
_phase_interpolator : Callable[[float], ArrayLike]
_image_position_interpolator : Callable[[float], ArrayLike]
_amplification_factor_interpolator : Callable[[float], ArrayLike]

def __calculate_critical_value__():
    '''
    For a given scaling of the NFW profile, the function will calculate the critical value of the
    source position. It does this by checking over a range of source positions for where the last
    occurance of three images is.
    '''

    source_position_space = np.linspace(0, 3.0, 1000)

    for idx, source_position in enumerate(source_position_space):
        image_positions = find_image_positions(source_position)
        number_of_images = len(image_positions)

        if number_of_images == 1:
            global _CRITICAL_VALUE
            _CRITICAL_VALUE = source_position_space[idx-1]
            break

def __generate_phase_interpolator__():
    '''
    For a given scaling of the NFW profile, the function will generate an interpolator that will
    return the phase required for the minimum time delay to be zero for given source positions.
    This interpolator is valid in the source position space below the critical value of the source
    psoition.
    '''

    source_position_space = np.linspace(0, _CRITICAL_VALUE, 1000)
    phase_space = np.empty(len(source_position_space))

    for idx, source_position in enumerate(source_position_space):
        phase_space[idx] = min_time_delay_phase(source_position)

    global _phase_interpolator
    _phase_interpolator = interp1d(source_position_space, phase_space)

def __generate_image_position_interpolator__():
    '''
    For a given scaling of the NFW profile, the function will generate an inteprolator that will
    return the image positions for the given value of source position. This interpolator is valid
    in the space below the critical value where three such positions will be generated.
    '''

    source_position_space = np.linspace(0, _CRITICAL_VALUE, 1000)
    image_one = np.empty(len(source_position_space))
    image_two = np.empty(len(source_position_space))
    image_three = np.empty(len(source_position_space))

    for idx, source_position in enumerate(source_position_space):
        image_positions = find_image_positions(source_position)
        image_one[idx] = image_positions[0]
        image_two[idx] = image_positions[1]
        image_three[idx] = image_positions[2]

    image_one_interpolator = interp1d(source_position_space, image_one)
    image_two_interpolator = interp1d(source_position_space, image_two)
    image_three_interpolator = interp1d(source_position_space, image_three)

    def image_position_interpolator(source_position):
        position_list = [image_one_interpolator(source_position),
                         image_two_interpolator(source_position),
                         image_three_interpolator(source_position)]
        return np.array(position_list)

    global _image_position_interpolator
    _image_position_interpolator = image_position_interpolator

def __generate_amplification_factor_interpolator__():
    '''
    For a given scaling of the NFW profile, the functin will generate an interpolator that will
    return for the amplification factor for the profile in the geometric optics regime.
    This interpolator is valid in the source position space above the critical value where only
    a single image is generated and the amplification is a constant value for that source position
    '''

    source_position_space = np.linspace(_CRITICAL_VALUE, 3.0, 1000)
    fiducial_dimensionless_frequency = 1000
    amplification_space = np.empty(len(source_position_space), dtype=complex)

    for idx, source_position in enumerate(source_position_space):
        amplification_space[idx] =\
            complete_amplification_factor(fiducial_dimensionless_frequency,
                                          source_position)

    amplification_real_space = np.real(amplification_space)
    amplification_imag_space = np.imag(amplification_space)

    real_interpolator = interp1d(source_position_space, amplification_real_space)
    imag_interpolator = interp1d(source_position_space, amplification_imag_space)

    def amplification_factor_interpolator(source_position):
        return real_interpolator(source_position) + 1j*imag_interpolator(source_position)

    global _amplification_factor_interpolator
    _amplification_factor_interpolator = amplification_factor_interpolator

def set_scaling(scaling_constant : float):
    '''
    Input:
        scaling_constant - characteristic scale length for the NFW profile

    For the specified scaling constant, the function sets the value of _SCALING and runs each
    of the subsequent functions that rely on the scaling constant value
    '''

    global _SCALING

    gravelogger.info("Setting NFW scaling value to %s", scaling_constant)
    _SCALING = scaling_constant

    __calculate_critical_value__()
    gravelogger.info("NFW Critical Value calulcated as %s", _CRITICAL_VALUE)

    __generate_phase_interpolator__()
    gravelogger.info("NFW Phase Interpolator Generated")

    __generate_image_position_interpolator__()
    gravelogger.info("NFW Image Position Interpolator Generated")

    __generate_amplification_factor_interpolator__()
    gravelogger.info("NFW Above Critical Amplification Factor Interpolator Generated")

def find_image_positions(source_position):
    '''
    Input:
        source_position - dimensionless displacement from the optical axis

    Output:
        position_array - array containing the positions of the images resulting from the lens
                         equation

    For the given source position and scaling constant, the function uses the C++ function
    ImagePositionArray to solve the lens equation and yield the positions of the lensed image. It
    then converts that to a numpy array for ease of python use.
    '''

    c_position_array = _cdll.PyImagePositions(ctypes.c_double(source_position),
                                              ctypes.c_double(_SCALING))

    position_array = np.ctypeslib.as_array(c_position_array, shape=(int(c_position_array[0]),))
    position_array = position_array[1:]

    print(position_array)

    return position_array

def min_time_delay_phase(source_position):
    '''
    Input:
        source_position - dimensionless displacement from the optical axis
        scaling_constant - characteristic scale length for the NFW profile

    Output:
        phase - value of the phase required for a minimum time delay of zero

    Function call the C++ backend function MinTimeDelayPhaseReal which calculates the value of the
    phase required for a minimum time delay of zero for the NFW profile at the given scale at the
    given source position.
    '''

    phase = float(_cdll.PyPhase(ctypes.c_double(source_position),
                                ctypes.c_double(_SCALING)))

    return phase

def complete_amplification_factor(dimensionless_frequency,
                                  source_position):
    '''
    Input:
        dimensionless_frequency - dimensionless form of the frequency of the gravitational wave data
        source_position - dimensionless displacement from the optical axis

    Output:
        amplification_value - complex value of the amplification factor

    Function uses the C++ backend function AFGRealOnly to compute the amplification factor for the
    NFW profile with the given input in the geometric optics regime. This version of the function
    requires no additional information, instead calculating the image positions and phase values
    directly in the function which makes it more complete but will make it more computationally
    expensive than the simpler version contained within SimpleAmpFac.
    '''

    c_result = _cdll.PyAmplificationFactorGeometric(ctypes.c_double(dimensionless_frequency),
                                                    ctypes.c_double(source_position),
                                                    ctypes.c_double(_SCALING),
                                                    None,
                                                    ctypes.c_double(0.0),
                                                    ctypes.c_int(0))
    amplification_value = complex(c_result[0], c_result[1])

    _cdll.destroyObj(c_result)

    return amplification_value

def amplification_factor(dimensionless_frequency_array,
                         source_position):
    '''
    Input:
        dimensionless_frequency_array - array containing the dimensionless form of the frequency
                                        of the gravitational wave data
        source_position - dimensionless displacement from the optical axis

    Output:
        amplification_array - complex values of the amplification factor for each dimensionless
                              frequency

    Function calculates the value of the amplification factor for the NFW profile in the geometric
    optics regime. To make this as efficient as possible, if the source position is above
    _CRITICAL_VALUE it will use the amplification factor interpolator that was precomputed when
    the scaling was set. If it is below that value, it will get the image positions and phase
    from the pregenerated interpolators to then feed into the C++ function SimpleAmpFac which
    allows for speedier computation of the value than the complete method outlined above.
    '''

    amplification_array = np.empty(len(dimensionless_frequency_array), dtype=complex)

    if source_position > _CRITICAL_VALUE:
        amplification_array[:] = _amplification_factor_interpolator(source_position)
        return amplification_array

    image_positions = _image_position_interpolator(source_position)
    number_of_images = len(image_positions)
    phase = _phase_interpolator(source_position)

    for idx, dimensionless_frequency in enumerate(dimensionless_frequency_array):
        c_result = _cdll.PyAmplificationFactorGeometric(
            ctypes.c_double(dimensionless_frequency),
            ctypes.c_double(source_position),
            ctypes.c_double(_SCALING),
            image_positions.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.c_double(phase),
            ctypes.c_int(number_of_images))

        amplification_array[idx] = complex(c_result[0], c_result[1])
        _cdll.destroyObj(c_result)

    return amplification_array

def generate_interpolator_data(config, args, file_dict):
    '''
    Input:
        config - INI configuration parser
        args - Commandline arguments passed to the program
        file_dict - dictionary of files containing dimensionless frequency and source poisition
                    values over which to generate the interpolator, followed by the corresponding
                    files containing the real and imaginary amplification factor values to use as
                    the interpolating data

    Function uses the C++ backend function within libnfw to generate the amplification factor data
    files from the input dimensionless frequency and source position files
    '''

    additional_arguments = get_additional_arguments(config, args,
                                                    _additional_arguments,
                                                    _additional_argument_types)

    gravelogger.info("Generating Lens Interpolator Data")
    _cdll.GenerateLensData(ctypes.c_char_p(file_dict["dimensionless_frequency"].encode("utf-8")),
                           ctypes.c_char_p(file_dict["source_position"].encode("utf-8")),
                           ctypes.c_char_p(file_dict["amplification_factor_real"].encode("utf-8")),
                           ctypes.c_char_p(file_dict["amplification_factor_imag"].encode("utf-8")),
                           ctypes.c_double(additional_arguments[0]),
                           ctypes.c_double(additional_arguments[1]),
                           ctypes.c_int(additional_arguments[2]),
                           ctypes.c_int(additional_arguments[3]))
    gravelogger.info("Lens Interpolator Data Generated")