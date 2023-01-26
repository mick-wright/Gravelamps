"""Model Interpolator Lensing Functions

Following are functions that perform calculations using a lensing interpolator generated from user
specified files. Included are the necessary functions to run the generators for such files for a
given model assuming it is similar in form to those used within Gravelamps directly.

Written by Mick Wright 2022

Globals
---------
_complex_interpolator
    Store of the generated complex interpolator
_lens_parameters
    Parameters to be passed to data generation functions

Routines
--------
amplification_factor
    Calculates amplification factor value for given grid values
generate_interpolator
    Generates interpolator from the specified files
"""

from collections.abc import Callable
import warnings

import numpy as np
from numpy.typing import ArrayLike
import scipy.interpolate as scint

_complex_interpolator : Callable[[float, float], ArrayLike]
_lens_parameters = ["lens_mass", "lens_fractional_distance", "source_position"]

def amplification_factor(dimensionless_frequency_array, source_position):
    """
    Calculates amplification factor value for given grid values using the generated interpolator.

    Parameters
    ----------
    dimensionless_frequency_array : float or Array of floats
        Dimensionless form of the frequency of gravitational wave data
    source_position : float
        Dimensionless displacement from the optical axis

    Returns
    -------
    amplification_array : float or Array or floats
        Interpolated values of the amplification factor for the specified grid values

    Raises
    ------
    AttributeError
        In case where the amplification factor is called before the interpolator is generated
    """

    if "_complex_interpolator" not in globals():
        raise AttributeError(
            "Interpolator has not been constructed, make sure generate_interpolator is run")

    amplification_array = _complex_interpolator(dimensionless_frequency_array, source_position)

    return amplification_array

def generate_interpolator(dimensionless_frequency_file,
                          source_position_file,
                          amplification_factor_real_file,
                          amplification_factor_imag_file):
    """
    Generates the amplification factor interpolator from the specified files.

    Specifically 1d interpolators are constructed for the real and imaginary components and these
    are combined into a single complex interpolator from the files specified. The result is placed
    into the _complex_interpolator global

    Parameters
    ----------
    dimensionless_frequency_file : str
        Path to file containing dimensionless frequency values to form grid axis
    source_position_file : str
        Path to file containing source position values to form grid axis
    amplification_factor_real_file : str
        Path to file containing real amplification factor data to form the interpolator
    amplification_factor_imag_file : str
        Path to file containing imaginary amplification factor data to form the interpolator
    """

    #Load in th files, creating arrays containing the data
    dimensionless_frequency_array = np.loadtxt(dimensionless_frequency_file)
    source_position_array = np.loadtxt(source_position_file)
    amplification_factor_real_array = np.loadtxt(amplification_factor_real_file)
    amplification_factor_imag_array = np.loadtxt(amplification_factor_imag_file)

    #Check that the arrays are the correct orientation --- assuming non-square matrix.
    #Issues a warning if the matrix is square
    if len(source_position_array) == len(dimensionless_frequency_array):
        warnings.warn(
            "Arrays are of equal length, cannot determine if matrix is correctly oriented")
    else:
        if amplification_factor_real_array.shape ==\
        (len(source_position_array), len(dimensionless_frequency_array)):
            amplification_factor_real_array = np.transpose(amplification_factor_real_array)
        if amplification_factor_imag_array.shape ==\
        (len(source_position_array), len(dimensionless_frequency_array)):
            amplification_factor_imag_array = np.transpose(amplification_factor_imag_array)

    #Construct the real and imaginary interpolators
    real_interpolator = scint.RectBivariateSpline(
        dimensionless_frequency_array, source_position_array, amplification_factor_real_array,
        kx=1, ky=1, s=0)
    imag_interpolator = scint.RectBivariateSpline(
        dimensionless_frequency_array, source_position_array, amplification_factor_imag_array,
        kx=1, ky=1, s=0)

    #Construct full complex interpolator
    def complex_interpolator(dimensionless_frequency, source_position):
        real_value = real_interpolator(dimensionless_frequency, source_position)
        imag_value = imag_interpolator(dimensionless_frequency, source_position)
        complex_value = real_value + 1j*imag_value
        return complex_value.flatten()

    global _complex_interpolator
    _complex_interpolator = complex_interpolator
