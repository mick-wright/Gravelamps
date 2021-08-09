'''
Helper functions for gravelamps.inference.

Functions generally aid in the process of getting the inference runs set up.
These are computational not physical in the nature of the work they do

Written by Mick Wright 2021
'''

import subprocess

from configparser import ConfigParser

import numpy as np

import gravelamps.lensing

def wy_handler(config):
    '''
    Input:
        config - INI configuration parser

    Output:
        dim_freq_file - location of the file containing the dimensionless frequency values
        sour_pos_file - location of the file containing the source position values

    Function handles the files containing the dimensionless frequency and source position values
    over which the amplification factor values will be calculated. If the user has specified files
    it will read them in, or will call the generator functions to make them.
    '''

    #Get the data subdirectory location
    outdir = config.get("output_settings", "outdir")
    data_subdirectory = outdir + "/data"

    #If the user has given locations of files containing the data, get those
    dim_freq_file = config.get(
        "lens_generation_settings", "dimensionless_frequency_file", fallback="None")
    sour_pos_file = config.get("lens_generation_settings", "source_position_file", fallback="None")

    #If the files are not given, call the generator functions
    if sour_pos_file == "None":
        sour_pos_file = gravelamps.lensing.utils.generate_source_position_file(config)
    if dim_freq_file == "None":
        dim_freq_file = gravelamps.lensing.utils.generate_dimensionless_frequency_file(config)

    #If external files are given, and the user has specfied, copy the files to the data subdirectory
    if config.getboolean("lens_generation_settings", "copy_files_to_data_subdirectory") is True:
        subprocess.run(["cp", dim_freq_file, data_subdirectory+"/w.dat"], check=True)
        subprocess.run(["cp", sour_pos_file, data_subdirectory+"/y.dat"], check=True)

        dim_freq_file = data_subdirectory+"/w.dat"
        sour_pos_file = data_subdirectory+"/y.dat"

    return dim_freq_file, sour_pos_file

def get_additional_parameters(config, geometric_optics_switch=0):
    '''
    Input:
        config - INI configuration parser
        geometric_optics_switch - optional value determining where in the array to switch from the
                                  wave optics calculation to the geometric optics approximation

    Output:
        additional_parameter_list - list containing the additional parameters needed for the model

    Function gets the additional parameters necessary for the function call necessary for the model.
    '''

    #Get the lens model
    lens_model = config.get("lens_generation_settings", "lens_model")

    #Get the precision
    precision = config.get("lens_generation_settings", "arithmetic_precision")

    #For each of the lens models, acquire the necessary additional parameters
    if lens_model == "pointlens":
        additional_parameter_list = [precision, geometric_optics_switch]

    elif lens_model == "sislens":
        integration_upper_limit = config.get(
            "lens_generation_settings", "sis_integration_upper_limit")
        additional_parameter_list = [
            integration_upper_limit, precision, geometric_optics_switch]

    elif lens_model == "nfwlens":
        scaling_constant = config.get(
            "lens_generation_settings", "nfw_scaling_constant")
        integration_upper_limit = config.get(
            "lens_generation_settings", "nfw_integration_upper_limit")
        additional_parameter_list = [
            scaling_constant, integration_upper_limit, precision, geometric_optics_switch]

    # If using a non-standard program, user may specify the parameters manually. Or may add
    # additional elif statements in this function
    else:
        additional_parameters = config.get(
            "lens_generation_settings", "additional_parameter_list").split(",")
        additional_parameter_list = []

        for parameter in additional_parameter_list:
            additional_parameter_list.append(config.get(
                "lens_generation_settings", parameter))

    return additional_parameter_list


def amp_fac_handler(config, dim_freq_file, sour_pos_file, mode="local"):
    '''
    Input:
        config - INI configuration parser
        dim_freq_file - location of the file containing the dimensionless frequency values
        sour_pos_file - location of the file containing the source position values
        mode - "local" or "pipe" determines whether to locally generate the data or make the
                condor submit files to generate the amplification factor data

    Outut:
        amp_fac_real_file - location of the file containing the real part of the amplification
                            factor value
        amp_fac_imag_file - location of the file containing the imaginary part of the amplification
                            factor value

    Function handles the amplification factor data for the interpolator. If the user has specified
    files it will read them, or it will generate the data/create the condor submit file
    '''

    #Get the data subdirectory location
    outdir = config.get("output_settings", "outdir")
    data_subdirectory = outdir + "/data"

    #If the user has given the locations of files containing the data, get those
    amp_fac_complex_file = config.get(
        "lens_generation_settings", "amplification_factor_complex_file", fallback="None")
    amp_fac_real_file = config.get(
        "lens_generation_settings", "amplification_factor_real_file", fallback="None")
    amp_fac_imag_file = config.get(
        "lens_generation_settings", "amplification_factor_imag_file", fallback="None")

    #If the complex file exists, read it in and then split it into real and imaginary parts
    if amp_fac_complex_file != "None":
        complex_array = np.loadtxt(amp_fac_complex_file, dtype=complex)
        real_array = np.real(complex_array)
        imag_array = np.imag(complex_array)

        amp_fac_real_file = data_subdirectory + "/fReal.dat"
        amp_fac_imag_file = data_subdirectory + "/fImag.dat"

        np.savetxt(amp_fac_real_file, real_array)
        np.savetxt(amp_fac_imag_file, imag_array)

    #Else if both of the real and imaginary files exist, bypass remaining functions
    elif amp_fac_real_file != "None" and amp_fac_imag_file != "None":
        pass

    #If the files do not exist, generate the data or the condor submission file
    else:
        amp_fac_real_file = data_subdirectory + "/fReal.dat"
        amp_fac_imag_file = data_subdirectory + "/fImag.dat"

        #Get the Lens Model
        lens_model = config.get("lens_generation_settings", "lens_model")

        #For the geometric optics switch, create the value by loading in the dimensionles frequency
        #file and finding the first place in the array where the array is greater than the
        #specified changeover frequency. In the case where the area is entirely below the switch
        #value, set the switch to the end of the array
        geometric_optics_min_dim_freq = config.getfloat(
            "lens_generation_settings", "geometric_optics_changeover_dimensionless_frequency")
        dim_freq_array = np.loadtxt(dim_freq_file)
        switch_value = np.argmax(dim_freq_array >= geometric_optics_min_dim_freq)
        if dim_freq_array[switch_value] < geometric_optics_min_dim_freq:
            switch_value = len(dim_freq_array)
        switch_value = str(switch_value)

        #Get the additional parameters necessary for the process call
        additional_parameters = get_additional_parameters(config, switch_value)
