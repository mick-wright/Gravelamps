'''
Gravelamps

Clusterised machine lens data generation and analysis.
Assumes cluster is using the HTCondor Scheduler

Written by Mick Wright 2021
'''

import os
import sys
import subprocess

from configparser import ConfigParser

import gravelamps.lensing
import gravelamps.inference

def main():
    '''
    Main function - takes the user generated INI to compile condor submission and DAG files
    in order to generate the lens data and perform the analysis run
    '''

    #Instantiate the Configuration Parser
    config = ConfigParser()

    #If user hasn't given a useable INI, raise exception
    if not os.path.isfile(sys.argv[1]):
        raise IOError("Input INI file not given!")

    #Check that the Configuration Parser can read the INI
    try:
        config.read(sys.argv[1])
    except IOError:
        print("Input INI file could not be read!")

    #Create the Outdir Directory and the Data and Submit Subfolders
    outdir = config.get("output_settings", "outdir")
    data_subdirectory = outdir + "/data"
    submit_subdirectory = outdir + "/submit"

    for folder in (outdir, data_subdirectory, submit_subdirectory):
        if not os.path.isdir(folder):
            os.mkdir(folder)

    #Generate hte Lensing Data
    dim_freq_file, sour_pos_file = gravelamps.inference.helpers.wy_handler(config)

    #Generate the Amplification Factor Submit File
    amp_fac_real_file, amp_fac_imag_file = gravelamps.inference.helpers.amp_fac_handler(
        config, dim_freq_file, sour_pos_file, mode="pipe")

    #Construct the Waveform Arguments Dictionary
    waveform_arguments = dict()

    waveform_approximant = config.get("analysis_settings", "waveform_approximant")
    minimum_frequency = config.getfloat("analysis_settings", "minimum_frequency")
    maximum_frequency = config.getfloat("analysis_settings", "maximum_frequency")
    reference_frequency = config.getfloat("analysis_settings", "reference_frequency")

    waveform_arguments["waveform_approximant"] = waveform_approximant
    waveform_arguments["minimum_frequency"] = minimum_frequency
    waveform_arguments["maximum_frequency"] = maximum_frequency
    waveform_arguments["reference_frequency"] = reference_frequency
    waveform_arguments["dim_freq_file"] = dim_freq_file
    waveform_arguments["sour_pos_file"] = sour_pos_file
    waveform_arguments["amp_fac_real_file"] = amp_fac_real_file
    waveform_arguments["amp_fac_imag_file"] = amp_fac_imag_file

    #If Injecting, Generate the Injection File and the Injection Waveform Arguments
    if config.getboolean("injection_settings", "injection"):
        #Read in the Injection Parameters, converting to floats
        injection_parameters = config._sections["injection_parameters"]
        injection_parameters.update(
            (key,float(value)) for key, value in injection_parameters.items())

        inject_file = gravelamps.inference.file_generators.injection_file(
            config, injection_parameters)

        #Construct the Injection Waveform Arguments dictionary
        injection_waveform_arguments = waveform_arguments.copy()

        dim_freq_other = config.get("injection_settings", "dimensionless_frequency_file")
        sour_pos_other = config.get("injection_settings", "source_position_file")
        amp_fac_real_other = config.get("injection_settings", "amp_fac_real_file")
        amp_fac_imag_other = config.get("injection_settings", "amp_fac_imag_file")

        if dim_freq_other != None:
            injection_waveform_arguments["dim_freq_file"] = os.path.abspath(dim_freq_other)
        if sour_pos_other != None:
            injection_waveform_arguments["sour_pos_file"] = os.path.abspath(sour_pos_other)
        if amp_fac_real_other != None:
            injection_waveform_arguments["amp_fac_real_file"] = os.path.abspath(amp_fac_real_other)
        if amp_fac_imag_other != None:
            injection_waveform_arguments["amp_fac_imag_file"] = os.path.abspath(amp_fac_imag_other)
    else:
        inject_file = None
        injection_waveform_arguments = None

    #If user specifies, perform unlensed analysis run
    if config.getboolean("unlensed_analysis_settings", "unlensed_analysis_run"):
        #Generate the unlensed run INI
        unlensed_ini = gravelamps.inference.file_generators.bilby_pipe_ini(
            config=config, inject_file=inject_file,
            injection_waveform_arguments=injection_waveform_arguments,
            waveform_arguments=waveform_arguments, mode="unlensed")

        #Run the bilby_pipe initial run
        subprocess.run(["bilby_pipe", unlensed_ini], check=True)

        #TODO: MODIFY PRIORS BASED ON UNLENSED RUN

    #Generate Lensed Run INI
    lensed_ini = gravelamps.inference.file_generators.bilby_pipe_ini(
        config=config, inject_file=inject_file,
        injection_waveform_arguments=injection_waveform_arguments,
        waveform_arguments=waveform_arguments, mode="lensed")

    #Run the bilby_pipe intial run
    subprocess.run(["bilby_pipe", lensed_ini], check=True)

    #Generate the overarching DAG file
    dag_file = gravelamps.inference.file_generators.overarching_dag(config)

    #Message user with submission
    print("To submit, use \n $ condor_submit_dag " + dag_file)
