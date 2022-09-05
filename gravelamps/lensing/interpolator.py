'''
Model Agnostic Interpolator Lensing Functions

These functions perform calculations using a lensing interpolator generated from user specified
files. The module contains the necessary functions to run the generators for the model assuming
it is either a python path or an executable program.

Written by Mick Wright 2022
'''

from collections.abc import Callable
import warnings

import numpy as np
from numpy.typing import ArrayLike
import scipy.interpolate as scint

_complex_interpolator : Callable[[float, float], ArrayLike]
_lens_parameters = ["lens_mass", "lens_fractional_distance", "source_position"]

def amplification_factor(dimensionless_frequency_array, source_position):
    '''
    Input:
        dimensionless_frequency_array - array containing the dimensionless form of the frequency
                                        of the gravitational wave data
        source_position - dimensionless displacement from the optical axis

    Output:
        amplification_array - complex values of the amplification factor for each dimensionless
                              frequency

    Function calculates the value of the amplification factor using the interpolator generated by
    generate_interpolator.
    '''

    if "_complex_interpolator" not in globals():
        raise AttributeError(
            "Interpolator has not been constructed, make sure generate_interpolator is run")

    amplification_array = _complex_interpolator(dimensionless_frequency_array, source_position)

    return amplification_array

def generate_interpolator(dimensionless_frequency_file,
                          source_position_file,
                          amplification_factor_real_file,
                          amplification_factor_imag_file):
    '''
    Input:
        dimensionless_frequency_file - string containing path to file containing dimensionless
                                       frequency values over which to generate the interpolator
        source_position_file - string containing path to file containing source position values
                               over which to generate the inteprolator
        amplification_factor_real_file - string containing path to file containing real component
                                         of the amplification factor over which to generate the
                                         interpolator
        amplification_factor_imag_file - string containing path to file containing imaginary
                                         component of the amplification factor over which to
                                         generate the interpolator

    Function takes in files containing an array of dimensionless frequency and source position as
    well as the corresponding amplification factor data which is used to generate interpolators
    for the real and imaginary components of the data which is combined into the single final
    complex interpolator returned.
    '''

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
