'''
Helper functions for gravelamps.inference.

Functions generally aid in the process of getting the inference runs set up.
These are computational not physical in the nature of the work they do

Written by Mick Wright 2021
'''

import os
import subprocess
import importlib

import numpy as np
import bilby

import gravelamps.lensing
import gravelamps.inference

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
        if dim_freq_file != data_subdirectory + "/w.dat":
            subprocess.run(["cp", dim_freq_file, data_subdirectory+"/w.dat"], check=True)
            dim_freq_file = data_subdirectory+"/w.dat"
        if sour_pos_file != data_subdirectory + "/y.dat":
            subprocess.run(["cp", sour_pos_file, data_subdirectory+"/y.dat"], check=True)
            sour_pos_file = data_subdirectory + "/y.dat"

    dim_freq_file = os.path.abspath(dim_freq_file)
    sour_pos_file = os.path.abspath(sour_pos_file)

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
            "lens_generation_settings", "sis_summation_upper_limit")
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
            "lens_generation_settings", "additional_parameter_list").replace(" ","").split(",")
        additional_parameter_list = []

        for parameter in additional_parameters:
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

        #Either directly call the lensing function if the mode is local or generate the submit file
        #if the mode is pipe
        if mode == "local":
            bilby.core.utils.logger.info("Generating Lens Data")
            proc_to_run = [lens_model, dim_freq_file, sour_pos_file, amp_fac_real_file,
                           amp_fac_imag_file] + additional_parameters
            subprocess.run(proc_to_run, check=True)
            bilby.core.utils.logger.info("Lens Data Generated")
        elif mode == "pipe":
            bilby.core.utils.logger.info("Generating Lens Data Submit File")
            gravelamps.inference.file_generators.lens_subfile(
                config, dim_freq_file, sour_pos_file, amp_fac_real_file,
                amp_fac_imag_file, additional_parameters)
            bilby.core.utils.logger.info("Lens Data Submit File Generated")

    # If the user has specified that the files be copied to the data subdirectory, do this
    if config.getboolean("lens_generation_settings", "copy_files_to_data_subdirectory"):
        if amp_fac_real_file != data_subdirectory + "/fReal.dat":
            subprocess.run(["cp", amp_fac_real_file, data_subdirectory+"/fReal.dat"], check=True)
            amp_fac_real_file = data_subdirectory + "/fReal.dat"
        if amp_fac_imag_file != data_subdirectory + "/fImag.dat":
            subprocess.run(["cp", amp_fac_imag_file, data_subdirectory+"/fImag.dat"], check=True)
            amp_fac_imag_file = data_subdirectory + "/fImag.dat"

    amp_fac_real_file = os.path.abspath(amp_fac_real_file)
    amp_fac_imag_file = os.path.abspath(amp_fac_imag_file)

    return amp_fac_real_file, amp_fac_imag_file

def wfgen_fd_source(waveform_generator_class_name, frequency_domain_source_model_name):
    '''
    Input:
        waveform_generator_class_name - string containing either the name of a bilby waveform
                                        generator class, or the full python path of an arbitrary
                                        waveform generator class
        frequency_domain_source_model_name - string containing either the name of a bilby frequency
                                             domain source model or the full pythoin path of an
                                             arbitrary frequency domain source model function

    Output:
        waveform_generator_class - the class specified by the input class name
        frequency_domain_source_model - the function specified by the input model name

    Function takes in strings containing the names of the waveform generator class and the
    frequency domain source model functions, checks if they're usable and returns the class and
    function to the user. Based upon the implementation in bilby_pipe
    '''

    #Waveform Generator - if it's a single name, assume part of the standard bilby class
    #Otherwise split the name to get the module and class names
    if "." in waveform_generator_class_name:
        waveform_generator_split = waveform_generator_class_name.split(".")
        waveform_generator_module = ".".join(waveform_generator_split[:-1])
        waveform_generator_class_name = waveform_generator_split[-1]
    elif waveform_generator_class_name in ("default", ""):
        waveform_generator_module = "bilby.gw.waveform_generator"
        waveform_generator_class_name = "WaveformGenerator"
    else:
        waveform_generator_module = "bilby.gw.waveform_generator"

    #Attempt to find class
    try:
        waveform_generator_class = getattr(
            importlib.import_module(
                waveform_generator_module), waveform_generator_class_name)
    except ImportError:
        print("Waveform Generator Class could not be found!")

    #Frequency Domain Source Model - if it's a single name, check if it's in the standard
    #bilby list, otherwise follow the same procedure as for the waveform generator
    if frequency_domain_source_model_name in bilby.gw.source.__dict__.keys():
        frequency_domain_source_model = bilby.gw.source.__dict__[
            frequency_domain_source_model_name]
    else:
        frequency_domain_source_model_split = frequency_domain_source_model_name.split(".")
        frequency_domain_source_model_module = ".".join(frequency_domain_source_model_split[:-1])
        frequency_domain_source_model_name = frequency_domain_source_model_split[-1]

        try:
            frequency_domain_source_model = getattr(
                importlib.import_module(
                    frequency_domain_source_model_module), frequency_domain_source_model_name)
        except ImportError:
            print("Frequency Domain Source Model Function could not be loaded!")

    return waveform_generator_class, frequency_domain_source_model
