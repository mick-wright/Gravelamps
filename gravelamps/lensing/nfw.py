"""Navarro, Frenk, White (NFW) Lensing Functions

Following are functions performing calculations for a Navarro, Frenk, White lens mass density
profile model. The module backend is based in libnfw.

Written by Mick Wright 2022

Globals
-------
_cdll : ctypes.CDLL
    Library of C++ functions necessary for calculations
_additional_arguments : list of str
    Additional arguments required to construct interpolator data
_additional_argument_types : list of types
    Types of the arguments that are given above
_lens_parameters : list of str
    Parameters used for the model

_SCALING : float
    Scaling constant used in the NFW profile
_CRITICAL_VALUE : float
    Critical value of source position for scale.

_phase_interpolator : Callable[[float], ArrayLike]
    Interpolating function to calculate morse phase for source position
_image_position_interpolator : Callable[[float], ArrayLike]
    Interpolating function to calculate the image position for source position
_amplification_factor_interpolator : Callable[[float], ArrayLike]
    Interpolating function for the region of amplification factor that is a single value i.e. above
    the critical value of source position

Routines
--------
__calculate_critical_value__
    Calculates critical value of source position
__generate_phase_interpolator__
    Generate interpolator for phase required for zero minimum time delay
__genrate_image_position_interpolator__
    Generate interpolator for the image value. Valid only below the critical value
__generate_amplification_factor_interpolator__
    Generate interpolator for the amplification factor in the single value regime. Valid only
    above the critical value of the source position.
set_scaling
    Set the value of _SCALING and regenerate the critical value and the interpolator functions
find_image_positions
    Calculates the position of the lensed images
min_time_delay_phase
    Calculats the value of the phase required for zero minimum time delay
complete_amplification_factor
    Calculates geometric optics amplification factor with no additional information
amplification_factor
    Calculates geometric optics amplification factor using the pregenerated interpolators where
    possible
generate_interpolator_data
    Generates the amplification factor data files for use in interpolator generation

"""

from collections.abc import Callable

import ctypes
import os

from importlib.resources import files

import numpy as np
from numpy.typing import ArrayLike
from scipy.interpolate import interp1d

from gravelamps.core.gravelog import gravelogger
from gravelamps.lensing.generic import get_additional_arguments

#The following loads in the DLL containing the C++ functions for the direct implementations for
#goemetric optics runs. It then sets the argument and result types accordingly

_cdll = ctypes.CDLL(importlib.resources.files('gravelamps.model.lib').joinpath('libnfw.so'))

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
    """
    Calculates critical value of source position.

    This calculation is done for the scaling of the NFW profile given in _SCALING. It is done by
    checking over a range of source positions for where the last occurance of three images is.
    occurance of three images is. It places the resultant value in _CRITICAL_VALUE.
    """

    source_position_space = np.linspace(0, 3.0, 1000)

    for idx, source_position in enumerate(source_position_space):
        image_positions = find_image_positions(source_position)
        number_of_images = len(image_positions)

        if number_of_images == 1:
            global _CRITICAL_VALUE
            _CRITICAL_VALUE = source_position_space[idx-1]
            break

def __generate_phase_interpolator__():
    """
    Generate interpolator for phase required for zero minimum time delay

    This calculation is done for the scaling of the NFW profile given in _SCALING. The interpolator
    generated is valid in the source position space below the critival value of the source position
    given in _CRITICAL_VALUE. The produced function is placed in _phase_interpolator.
    """

    source_position_space = np.linspace(0, _CRITICAL_VALUE, 1000)
    phase_space = np.empty(len(source_position_space))

    for idx, source_position in enumerate(source_position_space):
        phase_space[idx] = min_time_delay_phase(source_position)

    global _phase_interpolator
    _phase_interpolator = interp1d(source_position_space, phase_space)

def __generate_image_position_interpolator__():
    """
    Generate interpolator for the image value. Valid only below the critical value.

    This calculation is done for the scaling of the NFW profile given in _SCALING. The interpolator
    returns the three image positions for the given source position, making this interpolator valid
    only in the source position space below the critical value defined in _CRITICAL_VALUE where
    three images will be produced by the lensing. The produced function is placed in
    `_image_position_interpolator`.

    """

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
    """
    Generate interpolator for the amplification factor in the single value regime. Valid only above
    the critical value of source position.

    This calculation is done for the scaling of the NFW profile given in _SCALING. The interpolator
    will return a single value of amplification factor for any dimensionless frequency for a given
    source position i.e. this is valid in the regime where the amplification factor is flat. This
    is the case only above the critical value of source position defined in _CRITICAL_VALUE. The
    produced function is placed in _amplification_factor_interpolator
    """

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
    """
    Set the value of _SCALING and regenerate the critical value and interpolator functions.

    Parameters
    ----------
    scaling_constant : float
        Value to set the scaling constant to
    """

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
    """
    Calculates the position of the lensed images

    This calculation is done using the C++ function ImagePositionArray within libnfw. That function
    solves the lens equation to yield the positions of the lensed images.

    Parameters
    ----------
    source_position : float
        Dimensionless displacement from the optical axis

    Returns
    -------
    position_array : Array of floats
        Positions of the images resulting from the lens equation
    """

    c_position_array = _cdll.PyImagePositions(ctypes.c_double(source_position),
                                              ctypes.c_double(_SCALING))

    position_array = np.ctypeslib.as_array(c_position_array, shape=(int(c_position_array[0]),))
    position_array = position_array[1:]

    print(position_array)

    return position_array

def min_time_delay_phase(source_position):
    """
    Calculates the value of the phase required for zero minimum time delay

    This calculation is done using the C++ function MinTimeDelayPhaseReal within libnfw which
    calculates the real value of the phase required for a minimum time delay of zero for the given
    scaling.

    Parameters
    ----------
    source_position : float
        Dimensionless displacement from the optical axis

    Returns
    -------
    phase : float
        Phase required for zero minimum time delay
    """

    phase = float(_cdll.PyPhase(ctypes.c_double(source_position),
                                ctypes.c_double(_SCALING)))

    return phase

def complete_amplification_factor(dimensionless_frequency,
                                  source_position):
    """
    Calculates geometric optics amplification factor with no additional information

    The calculation is done using the C++ function AFGRealOnly within libnfw. This calculation
    does not require additional information beyond of a value of _SCALING. As a consequence it is
    always complete, but is not as computationally efficient compared to `amplification_factor`.

    Parameters
    ----------
    dimensionless_frquency : float
        Dimensionless form of the frequency of interest
    source_position : float
        Dimensionless displacement from the optical axis

    Returns
    -------
    amplification_value : complex
        Amplification of signal produced by lensing

    """

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
    """
    Calculates geometric optics amplification factor using the pregenerated information where
    possible.

    This calculation has been designed to be more efficient than `complete_amplification_factor` by
    using the pregnenerated interpolating functions. Above a source position of `_CRITICAL_VALUE`
    it will use the single value `_amplification_factor_interpolator`. Below, it will use the C++
    function `SimpleAmpFac` which will take in the image positions and phase computed from the
    relevant interpolators.

    Parameters
    ----------
    dimensionless_frquency_array : Array of floats
        Dimensionless form of the frequencies of interest
    source_position : float
        Dimensionless displacement from the optical axis

    Returns
    -------
    amplification_array : Array of complex
        Amplification of signal produced by lensing over specified grid
    """

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
    """
    Generates the amplification factor data files for use in interpolator generation

    This is done via the C++ `GenerateLensData` function within `libnfw`. It will read in the
    specified grid files and fill the data files with the appropraite values of amplification
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
                           ctypes.c_double(additional_arguments[0]),
                           ctypes.c_double(additional_arguments[1]),
                           ctypes.c_int(additional_arguments[2]),
                           ctypes.c_int(additional_arguments[3]))
    gravelogger.info("Lens Interpolator Data Generated")
