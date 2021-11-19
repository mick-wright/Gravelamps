'''
Utility functions for gravelamps.lensing

These functions perform the construction and generation of the lensing code
more directly.

Written by Mick Wright 2021
'''

import scipy.interpolate as scint
import numpy as np

def generate_dimensionless_frequency_file(config):
    '''
    Input:
        config - INI configuration parser

    Output:
        dim_freq_file - location of the generated dimensionless frequency file

    Function from the options set in the INI, generates a file called "w.dat" in which there will
    be an appropriately constructed list of dimensionless frequency values over which the
    interpolator will be generated
    '''

    # Construct the full path to the file
    outdir = config.get("output_settings", "outdir")
    dim_freq_file = outdir + "/data/w.dat"

    # Get settings from the INI file determining the minimum and maximum values to consider
    # as well as the number of points to generate
    min_value = config.getfloat("lens_generation_settings", "minimum_dimensionless_frequency")
    max_value = config.getfloat("lens_generation_settings", "maximum_dimensionless_frequency")
    num_values = config.getint("lens_generation_settings", "dimensionless_frequency_length")

    # Generate the dimensionless frequency array
    dim_freq_array = np.linspace(min_value, max_value, num_values)

    # Save the resultant array to the filename
    np.savetxt(dim_freq_file, dim_freq_array)

    return dim_freq_file

def generate_source_position_file(config):
    '''
    Input:
        config - INI configuration parser

    Output:
        sour_pos_file - location of the generated source position file

    Function from the options set in the INI, generates a file called "y.dat" in which there will
    be an appropriately constructed list of source position values over which the interpolaotr
    will be generated
    '''

    # Construct the full path to the file
    outdir = config.get("output_settings", "outdir")
    sour_pos_file = outdir + "/data/y.dat"

    # Get settings from the INI file determining the minimum and maximum values to consider
    # as well as the number of points to generate
    min_value = config.getfloat("lens_generation_settings", "minimum_source_position")
    max_value = config.getfloat("lens_generation_settings", "maximum_source_position")
    num_values = config.getint("lens_generation_settings", "source_position_length")

    # Generate the source position array
    sour_pos_array = np.linspace(min_value, max_value, num_values)

    # Save the resultant array to the filename
    np.savetxt(sour_pos_file, sour_pos_array)

    return sour_pos_file

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
