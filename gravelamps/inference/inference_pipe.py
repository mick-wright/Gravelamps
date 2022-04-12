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
    data_subdirectory = f"{outdir}/data"
    submit_subdirectory = f"{outdir}/submit"

    for folder in (outdir, data_subdirectory, submit_subdirectory):
        if not os.path.isdir(folder):
            os.mkdir(folder)

    #Get which methodology is being used
    methodology = config.get("lens_generation_settings", "methodology")

    #Dependent upon methodology selected, construct waveform arguments dictionary
    if methodology == "interpolate":
        #Generate the Lensing Data
        dim_freq_file, sour_pos_file = gravelamps.inference.helpers.wy_handler(config)

        #Generate the Amplification Factor Submit File
        amp_fac_real_file, amp_fac_imag_file = gravelamps.inference.helpers.amp_fac_handler(
            config, dim_freq_file, sour_pos_file, mode="pipe")

        #Actually construct the dictionary
        waveform_arguments = gravelamps.inference.helpers.construct_waveform_arguments(
            config, "analysis", dim_freq_file=dim_freq_file, sour_pos_file=sour_pos_file,
            amp_fac_real_file=amp_fac_real_file, amp_fac_imag_file=amp_fac_imag_file)

    elif methodology == "direct":
        waveform_arguments = gravelamps.inference.helpers.construct_waveform_arguments(
            config, "analysis")

    #If Injecting, Generate the Injection File and the Injection Waveform Arguments
    if config.getboolean("injection_settings", "injection"):
        #Read in the Injection Parameters, converting to floats
        injection_parameters = {}

        for key, value in config.items("injection_parameters"):
            injection_parameters[key] = float(value)

        inject_file = gravelamps.inference.file_generators.injection_file(
            config, injection_parameters)

        #Load in the injection methodology to see if it isn't None
        injection_methodology = config.get("injection_settings", "methodology")

        #Construct injection waveform arguments
        if injection_methodology == "None":
            injection_waveform_arguments = waveform_arguments.copy()
        else:
            injection_waveform_arguments =\
                gravelamps.inference.helpers.construct_waveform_arguments(config, "data")

    else:
        inject_file = None
        injection_waveform_arguments = None

    #Generate Lensed Run INI
    bilby_pipe_ini = gravelamps.inference.file_generators.bilby_pipe_ini(
        config=config, inject_file=inject_file,
        injection_waveform_arguments=injection_waveform_arguments,
        waveform_arguments=waveform_arguments)

    #Run the bilby_pipe intial run
    subprocess.run(["bilby_pipe", bilby_pipe_ini], check=True)

    #Generate the overarching DAG file
    dag_file = gravelamps.inference.file_generators.overarching_dag(config)

    #Message user with submission
    print(f"To submit, use \n $ condor_submit_dag {dag_file}")
