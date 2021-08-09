'''
Utility functions for gravelamps.lensing

These functions perform the construction and generation of the lensing code
more directly.

Written by Mick Wright 2021
'''

import numpy as np

from configparser import ConfigParser

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
    min_value = config.getint("lens_generation_settings", "minimum_dimensionless_frequency")
    max_value = config.getint("lens_generation_settings", "maximum_dimensionless_frequency")
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
    min_value = config.getint("lens_generation_settings", "minimum_source_position")
    max_value = config.getint("lens_generation_settings", "maximum_source_position")
    num_values = config.getint("lens_generation_settings", "source_position_length")

    # Generate the source position array
    sour_pos_array = np.linspace(min_value, max_value, num_values)

    # Save the resultant array to the filename
    np.savetxt(sour_pos_file, sour_pos_array)

    return sour_pos_file
