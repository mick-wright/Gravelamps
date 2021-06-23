'''
Utility functions for Gravelamps

Written by Mick Wright 2021
'''

import os
import subprocess
import importlib

from configparser import ConfigParser

import numpy as np
import scipy.interpolate as scint
import bilby

import gwlensing.lensing

def generate_dimensionless_frequency_file(config):
    '''
    Input:
        config - Ini configuartion parser

    Output:
        filename - containing the location of generated w file

    Function generates a file in the data subdirectory called "w.dat" which will
    contain a set of dimensionless frequency values over which the interpolation will occur
    '''

    #Get location for the data file to go
    outdir = config.get("bilby_setup", "outdir")
    data_subdir = outdir + "/" + config.get("data_settings", "data_subdir")
    filename = data_subdir + "/w.dat"

    #Acquire Data from ini file
    min_dim_freq = config.getint("data_settings", "min_dimensionless_frequency_for_w_generation")
    max_dim_freq = config.getint("data_settings", "max_dimensionless_frequency_for_w_generation")
    num_points_to_gen = config.getint("data_settings", "number_of_dimensionless_frequency_points")

    #Generate the dimensionless frequency array
    dim_freq_array = np.linspace(min_dim_freq, max_dim_freq, num_points_to_gen)

    #Save resultant array to w.dat data file
    np.savetxt(filename, dim_freq_array)

    return filename

def generate_impact_parameter_file(config):
    '''
    Input:
        config - Ini configuration parser

    Output:
        filename - containing the location of generated y file

    Function generates a file in the data subdirectory called "y.dat" which will
    contain a set of dimensionless frequency values over which the interpolation will occur
    '''

    #Get location for the data file to go
    outdir = config.get("bilby_setup", "outdir")
    data_subdir = outdir + "/" + config.get("data_settings", "data_subdir")
    filename = data_subdir + "/y.dat"

    #Acquire Data from Ini file
    min_y = config.getfloat("data_settings", "min_impact_parameter")
    max_y = config.getfloat("data_settings", "max_impact_parameter")
    num_points_to_gen = config.getint("data_settings", "number_of_impact_parameter_points")

    #Generate the impact parameter array
    y_array = np.linspace(min_y, max_y, num_points_to_gen)

    #Save resultant array to y.dat data file
    np.savetxt(filename, y_array)

    return filename

def wy_handler(config):
    '''
    Input:
        config - Ini configuration parser

    Output:
        w_array_file - location of the file containing the dimensionless frequecy values
        y_array_file - location of the file containing the impact parameter values

    Function handles making sure that the dimensionless frequency and impact parameters are
    correctly handled given and will move externally given files to the data subdir if
    the user specifies
    '''

    #Get location for the data files
    outdir = config.get("bilby_setup", "outdir")
    data_subdir = outdir + "/" + config.get("data_settings", "data_subdir")

    #Get optionally input user files
    w_array_file = config.get("optional_input", "dimensionless_frequency_file", fallback="None")
    y_array_file = config.get("optional_input", "impact_parameter_file", fallback="None")

    #For both dimensionless frequency and impact parameter
    #if optional input has not been given check if the file already exists in the data subdir
    #otherwise generate a new one
    if w_array_file == "None":
        if not os.path.isfile(data_subdir+"/w.dat"):
            w_array_file = gwlensing.lensing.utils.generate_dimensionless_frequency_file(config)
        else:
            w_array_file = data_subdir + "/w.dat"

    if y_array_file == "None":
        if not os.path.isfile(data_subdir+"/y.dat"):
            y_array_file = gwlensing.lensing.utils.generate_impact_parameter_file(config)
        else:
            y_array_file = data_subdir + "/y.dat"

    #If external files are used and the copy to data settings is enabled copy files to data subdir
    if config.getboolean("optional_input", "copy_files_to_data") is True:
        subprocess.run(["cp", w_array_file, data_subdir+"/w.dat"], check=True)
        subprocess.run(["cp", y_array_file, data_subdir+"/y.dat"], check=True)

    return w_array_file, y_array_file

def get_additional_parameters(config):
    '''
    Input:
        config - Ini configuration parser

    Output:
        additional_parameter_list - list containing the additional parameters needed for the model

    Function gets the additional parameters necessary for the function call necessary for the model
    '''

    #Get the lens model
    lens_model = config.get("lens_settings", "lens_model")

    #For each of the lens models, acquire the necessary additional parameters
    if lens_model == "pointlens":
        additional_parameter_list = []
    elif lens_model == "sislens":
        additional_parameter_list = []
    elif lens_model == "nfwlens":
        ks_val = config.get("lens_settings", "nfw_ks_val")
        integration_upper_limit = config.get("lens_settings", "integration_upper_limit")
        precision = config.get("lens_settings", "arithemtic_precision")

        additional_parameter_list = [ks_val, integration_upper_limit, precision]

    #If user wishes to add additional models here, simply include what additional parameters are
    #necessary here, naming the lens_model the same as the name of the function that does the
    #calculation

    else:
        additional_parameter_list = []

    return additional_parameter_list

def amp_fac_handler(config, w_array_file, y_array_file, mode="local"):
    '''
    Input:
        config - Ini configuration parser
        w_array_file - filename for the dimensionless frequency data for the interpolator
        y_array_file - filename for the impact parameter data for the interpolator
        mode - can be either "local" or "pipe" which is whether to run the command locally or
        generate a condor submit file for a cluster run respectively

    Output:
        amp_fac_real_file - filename for the real part of the generated amplification factor values
        amp_fac_imag_file - filename for the imaginary part of the generated amplification factor
        values

    Function handles the amplification factor data for the interpolator - either generating or
    reading it in and moving it depending upon user decisions.
    '''

    #Get the data subdirectory
    outdir = config.get("bilby_setup", "outdir")
    data_subdir = outdir + "/" + config.get("data_settings", "data_subdir")

    #Get the Lensing Model and the necessary additional parameters
    lens_model = config.get("lens_settings", "lens_model")
    additional_parameters = get_additional_parameters(config)

    #Read in the values of the optional input values
    amp_fac_complex_file = config.get("optional_input",
                                      "amplification_factor_complex_file", fallback="None")
    amp_fac_real_file = config.get("optional_input",
                                   "amplification_factor_real_file", fallback="None")
    amp_fac_imag_file = config.get("optional_input",
                                   "amplification_factor_imag_file", fallback="None")

    #If the complex file exists, read it in and split it into real and imaginary parts
    if amp_fac_complex_file != "None":
        complex_array = np.loadtxt(amp_fac_complex_file, dtype=complex)
        real_array = np.real(complex_array)
        imag_array = np.imag(complex_array)

        amp_fac_real_file = data_subdir + "/fReal.dat"
        amp_fac_imag_file = data_subdir + "/fImag.dat"

        np.savetxt(amp_fac_real_file, real_array)
        np.savetxt(amp_fac_imag_file, imag_array)

        #Else if both the real and imaginary files exists, read them in
    elif amp_fac_real_file != "None" and amp_fac_imag_file != "None":
        pass

        #Otherwise generate either the amplification files directly, or the submit file for condor
    else:
        amp_fac_real_file = data_subdir + "/fReal.dat"
        amp_fac_imag_file = data_subdir + "/fImag.dat"

        if mode == "local":
            print("Generating Lens Data")
            proc_to_run = [lens_model, w_array_file, y_array_file, amp_fac_real_file,
                           amp_fac_imag_file] + additional_parameters
            subprocess.run(proc_to_run, check=True)
            print("Lens Data Generated")
        elif mode == "pipe":
            print("Generating Lens Data Submit File")
            generate_lens_subfile(config, w_array_file, y_array_file,
                                  amp_fac_real_file, amp_fac_imag_file)
            print("Lens Data Submit File Generated")

    #If external files are used and the user has specified copy the files to the data subdirectory
    if config.getboolean("optional_input", "copy_files_to_data"):
        subprocess.run(["cp", amp_fac_real_file, data_subdir+"/fReal.dat"], check=True)
        subprocess.run(["cp", amp_fac_imag_file, data_subdir+"/fImag.dat"], check=True)

    return amp_fac_real_file, amp_fac_imag_file

def generate_interpolator(w_array_file, y_array_file, amp_fac_real_file, amp_fac_imag_file):
    '''
    Input:
        w_array_file - file containing the dimensionless frequency array over which to interpolate
        y_array_file - file containing the impact parameter array over which to interpolate
        amp_fac_real_file - file containing the real part of the amplification factor data
        amp_fac_imag_file - file containing the imaginary part of the amplification factor data

    Output:
        interpolator_func - a function which calculates the interpolation function

    Function takes in the dimensionless frequency and impact parameter arrays as well as the
    amplification factor data and constructs two interpolators for the real and imaginary parts
    and generates a final complex interpolator
    '''

    #From the files, load in the arrays
    w_array = np.loadtxt(w_array_file)
    y_array = np.loadtxt(y_array_file)
    amp_fac_real = np.loadtxt(amp_fac_real_file)
    amp_fac_imag = np.loadtxt(amp_fac_imag_file)

    #Make sure that the arrays are the correct orientation
    if amp_fac_real.shape == (len(y_array), len(w_array)):
        amp_fac_real = np.transpose(amp_fac_real)
    if amp_fac_imag.shape == (len(y_array), len(w_array)):
        amp_fac_imag = np.transpose(amp_fac_imag)

    #Construct the individual interpolators for the real and imaginary parts
    amp_interp_real = scint.RectBivariateSpline(w_array, y_array, amp_fac_real, kx=1, ky=1, s=0)
    amp_interp_imag = scint.RectBivariateSpline(w_array, y_array, amp_fac_imag, kx=1, ky=1, s=0)

    #Define the full complex interpolator for the amplification factor
    def interpolator_func(w_val, y_val):
        return (amp_interp_real(w_val, y_val) + 1j*amp_interp_imag(w_val, y_val)).flatten()

    return interpolator_func

def generate_lens_subfile(config, w_array_file, y_array_file, amp_fac_real_file, amp_fac_imag_file):
    '''
    Input:
        config - Ini configuration parser
        w_array_file - file containing the dimensionless frequency array
        y_array_file - file containing the impact parameter array over which to interpolate
        amp_fac_real_file - file containing the real part of the amplification factor
        amp_fac_imag_file - file contianing the imaginary part of the amplification factor data

    Function generates the condor submission file in order to generate the amplification
    factor data necessary for the full runs
    '''

    #Get the submission subdirectory
    outdir = config.get("bilby_setup", "outdir")
    submit_directory = outdir + "/submit"

    #Get the lens model
    lens_model = config.get("lens_settings", "lens_model")

    #Get the executable directory and construct the full path of the lens model function
    executable_directory = config.get("lens_settings", "executable_directory")
    executable_path = executable_directory + "/" + lens_model

    #Get the additional lens parameters needed for the lens models
    additional_lens_parameters = get_additional_parameters(config)

    #Create the submission directory if it doesn't exist
    if not os.path.isdir(submit_directory):
        os.mkdir(submit_directory)

    #Create the subfile and get the condor settings from the main ini file
    subfile = submit_directory + "/generate_lens.sub"
    condor_settings_dictionary = config._sections["condor_settings"].copy()

    #Open the submit file and write the necessary parts
    sub = open(subfile, "w")
    sub.write("universe = vanilla\n")
    sub.write("transfer_input_files = " + os.path.abspath(w_array_file) + ","
              + os.path.abspath(y_array_file) + "\n")
    sub.write("executable = " + executable_path + "\n")

    #Construct the arguments necessary for the function call
    arguments = ""

    for argument in (w_array_file, y_array_file, amp_fac_real_file, amp_fac_imag_file):
        to_add = argument + " "
        arguments += to_add
    for argument in additional_lens_parameters:
        to_add = argument + " "
        arguments += to_add

    sub.write("arguments = " + arguments + "\n")
    sub.write("log = " + submit_directory + "/lens_generation.log\n")
    sub.write("output = " + submit_directory + "/lens_generation.out\n")
    sub.write("error = " + submit_directory + "/lens_generation.err\n")

    for key, value in condor_settings_dictionary.items():
        sub.write(key + " = " + value + "\n")

    sub.write("queue 1\n")

def gen_inject_file(config, injection_parameters):
    '''
    Input:
        config - Ini configuration parser
        injection_parameters - dictionary containing the injection parameters

    Output:
        inject_dat_filename - string containing the path to the injection file

    Function generates the necessary injection file for performing analysis on simulated data
    '''

    print("Generating Injection File")

    #Get the proper location for the injection file
    outdir = config.get("bilby_setup", "outdir")
    inject_filename = outdir + "/data/injection.prior"
    inject_dat_filename = outdir + "/data/injection.dat"

    #Open the injection file, write the injection parameters to it, and then close the file
    inject_file = open(inject_filename, "w")

    for key, value in injection_parameters.items():
        inject_file.write(str(key) + "=" + str(value) + "\n")

    inject_file.close()

    #Run the bilby-pipe create_inject_file function
    n_injections = config.get("bilby_pipe_settings", "n-injections")

    subprocess.run(["bilby_pipe_create_injection_file", inject_filename, "--n-injection",
                    n_injections, "-f", inject_dat_filename], check=True)

    print("Injection File Generated")

    return inject_dat_filename

def gen_bilby_pipe_ini(config, inject_file, waveform_arguments, mode):
    '''
    Input:
        config - Ini configuration parser
        inject_file - Filename for the injection file
        waveform_arguments - Dictionary containing the waveform arguments
        mode - string can either be "lensed" or "unlensed" - tells the function which kind of run
               to generate

    Output:
        ini_filename - filename of the generated ini file

    Function creates a bilby_pipe style ini from the options given by the user in
    the gravelamps ini
    '''

    #Get the bilby_pipe ini filename and open it for writing
    bilby_pipe_ini_filename = config.get("bilby_setup", "label") + mode + "_bilby_pipe.ini"
    bilby_pipe_ini = open(bilby_pipe_ini_filename, "w")

    #Create the empty configuration dictionary
    bilby_pipe_config = dict()

    #Insert the accounting tag
    bilby_pipe_config["accounting"] = config.get("condor_settings", "accounting")

    #Insert the Label and the Output Directory
    bilby_pipe_config["label"] = config.get("bilby_setup", "label") + mode
    bilby_pipe_config["outdir"] = config.get("bilby_setup", "outdir")

    #Create the detector list
    detectors = config.get("bilby_setup", "detectors")
    detector_list = list(detectors.split(","))
    bilby_pipe_config["detectors"] = detector_list

    #Include the Waveform Generator Class and Frequency Domain Source Model
    if mode == "lensed":
        bilby_pipe_config["waveform-generator-class"] = config.get(
                "bilby_setup", "lensed_waveform_generator_class")
        bilby_pipe_config["frequency-domain-source-model"] = config.get(
                "bilby_setup", "lensed_frequency_domain_source_model")
    elif mode == "unlensed":
        bilby_pipe_config["waveform-generator-class"] = config.get(
                "data_settings", "unlensed_waveform_generator_class")
        bilby_pipe_config["frequency-domain-source-model"] = config.get(
                "data_settings", "unlensed_frequency_domain_source_model")

    #Include the duration and sampling frequency
    bilby_pipe_config["duration"] = config.get("bilby_setup", "duration")
    bilby_pipe_config["sampling_frequency"] = config.get("bilby_setup", "sampling_frequency")

    #Include settings from bilby_pipe_settings
    for key, value in config._sections["bilby_pipe_settings"].items():
        bilby_pipe_config[key] = value

    #If given, include injection file
    if inject_file is not None:
        bilby_pipe_config["injection_file"] = inject_file

    #Sampler Settings - Creating Dictionary from the settings
    bilby_pipe_config["sampler"] = config.get("bilby_setup", "sampler")
    bilby_pipe_config["sampler-kwargs"] = config._sections["sampler_settings"].copy()

    #Include the Prior File
    bilby_pipe_config["prior-file"] = config.get("prior_settings", "prior_file")

    #Include the waveform_arguments dictionary
    bilby_pipe_config["waveform_arguments_dict"] = waveform_arguments

    #Write the dictionary to the file
    for key, value in bilby_pipe_config.items():
        bilby_pipe_ini.write(key + " = " + value + "\n")

    #Close the ini file
    bilby_pipe_ini.close()

    return bilby_pipe_ini_filename

def gen_overarch_dag(config):
    '''
    Input:
        config - Ini configuration parser

    Function generates the overarching DAG file to submit all of the condor jobs -
    the lensing generation and the bilby pipe analysis runs
    '''

    #Get the submission directory
    outdir = config.get("bilby_setup", "outdir")
    submit_directory = outdir + "/submit"

    #Get the Lens Generation submit file
    lens_generation_subfile = submit_directory + "/generate_lens.sub"

    #Get the bilby-pipe dag file
    for filename in os.listdir(submit_directory):
        if filename.startswith("dag"):
            bilby_pipe_dag = submit_directory + "/" + filename

    #Get the filename for the overarching dag file
    label = config.get("bilby_setup", "label")
    overarch_filename = submit_directory + "/dag_" + label + "_overarch.submit"

    #Open the DAG file and write it
    overarch_dag = open(overarch_filename, "w")

    overarch_dag.write("JOB lens_generation " + os.path.abspath(lens_generation_subfile) + "\n")
    overarch_dag.write("SUBDAG EXTERNAL bilby_pipe " + os.path.abspath(bilby_pipe_dag) + "\n")
    overarch_dag.write("PARENT lens_generation CHILD bilby_pipe")

    overarch_dag.close()

def wfgen_fd_source(waveform_generator_class_name, frequency_domain_source_model_name):
    '''
    Input:
        waveform_generator_class_name - string containing either the name of a bilby
            waveform generator class or the full python path of an arbitrary
            waveform generator class
        frequency_domain_source_model_name - string containing either the name of a bilby
            frequency domain source model or the full python path of an arbitrary
            frequency domain source model function

    Output:
        waveform_generator_class - the class specified by the input class name
        frequency_domain_source_model - the function specified by the input model name

    Function takes in strings containing the names of the waveform generator class and the
    frequency domain source model function, checks if they're usable and returns the class
    and function to the user. Based upon the implementations in bilby-pipe
    '''

    #Waveform Generator - if it's a single name, assume part of standard bilby class
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

    #Now attempt to find the class
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
