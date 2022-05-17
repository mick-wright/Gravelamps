'''
Utility functions for gravelamps.lensing

These functions perform the construction and generation of the lensing code
more directly.

Written by Mick Wright 2021
'''

import scipy.interpolate as scint
import numpy as np

def generate_value_file(config, value_type, output_filename):
    '''
    Input:
        config - INI configuration parser
        value_type - Full name of the type of values to be genrated e.g dimensionless_frequency
        output_filename - filepath to the file to generate

    File generates the seuqence of values for the type value under consideration saving it to the
    specified filename
    '''

    min_value = config.getfloat("lens_interpolation_settings", f"minimum_{value_type}")
    max_value = config.getfloat("lens_interpolation_settings", f"maximum_{value_type}")
    num_values = config.getint("lens_interpolation_settings", f"{value_type}_length")

    value_array = np.linspace(min_value, max_value, num_values)

    np.savetxt(output_filename, value_array)

def generate_interpolator(dim_freq_file, sour_pos_file, amp_fac_real_file, amp_fac_imag_file):
    '''
    Input:
        dim_freq_file - location of file containing dimensionless frequency values to generate
                        interpolator over
        sour_pos_file - location of file containing source position values to generate interpolator
                        over
        amp_fac_real_file - location of file containing the real part of the amplification factor
                            values to interpolate over
        amp_fac_imag_file - location of file containing the imaginary part of the amplification
                            factor values to interpolate over

    Output:
        interpolator_func - function calculating the amplification factor at a given dimensionless
                            frequency and source position by interpolating over the input data

    Function takes in files containing an array of dimensionless frequency and source position as
    well as the amplification factor data for these values, constructs two interpolators for the
    real and imaginary parts of the value, and finally returns the complex interpolator giving the
    full value
    '''

    #From the files, load in the arrays
    dim_freq_array = np.loadtxt(dim_freq_file)
    sour_pos_array = np.loadtxt(sour_pos_file)
    amp_fac_real = np.loadtxt(amp_fac_real_file)
    amp_fac_imag = np.loadtxt(amp_fac_imag_file)

    #Make sure that the arrays are the correct orientation - assuming non square matrix
    if len(sour_pos_array) != len(dim_freq_array):
        if amp_fac_real.shape == (len(sour_pos_array), len(dim_freq_array)):
            amp_fac_real = np.transpose(amp_fac_real)
        if amp_fac_imag.shape == (len(sour_pos_array), len(dim_freq_array)):
            amp_fac_imag = np.transpose(amp_fac_imag)

    #Construct the real and imaginary interpolators
    real_interpolator = scint.RectBivariateSpline(
        dim_freq_array, sour_pos_array, amp_fac_real, kx=1, ky=1, s=0)
    imag_interpolator = scint.RectBivariateSpline(
        dim_freq_array, sour_pos_array, amp_fac_imag, kx=1, ky=1, s=0)

    #Construct the full complex interpolating function
    def interpolator_func(dim_freq, sour_pos):
        real_value = real_interpolator(dim_freq, sour_pos)
        imag_value = imag_interpolator(dim_freq, sour_pos)
        complex_value = real_value + 1j*imag_value
        return complex_value.flatten()

    return interpolator_func
