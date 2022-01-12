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

    #Get which methodology is being used
    methodology = config.get("analysis_settings", "methodology")

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

    elif methodology in ("direct-nonnfw", "direct-nfw"):
        waveform_arguments = gravelamps.inference.helpers.construct_waveform_arguments(
            config, "analysis", dim_freq_file=dim_freq_file, sour_pos_file=sour_pos_file,
            amp_fac_real_file=amp_fac_real_file, amp_fac_imag_file=amp_fac_imag_file)

    #If Injecting, Generate the Injection File and the Injection Waveform Arguments
    if config.getboolean("injection_settings", "injection"):
        #Read in the Injection Parameters, converting to floats
        injection_parameters = config._sections["injection_parameters"]
        injection_parameters.update(
            (key,float(value)) for key, value in injection_parameters.items())

        inject_file = gravelamps.inference.file_generators.injection_file(
            config, injection_parameters)

        #Load in the injection methodology to see if it isn't None
        injection_methodology = config.get("injection_settings", "methodology")

        #Construct injection waveform arguments
        if injection_methodology is None:
            injection_waveform_arguments = waveform_arguments.copy()
        else:
            injection_waveform_arguments = (
                gravelamps.inference.helpers.construct_waveform_arguments(config, "data"))

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
