'''
Gravelamps Pipe Inference

Condor based lens generation and analysis variant
for Simulated Data

Written by Mick Wright 2021
'''

import os
import sys
import subprocess

from configparser import ConfigParser

import bilby

import gwlensing.lensing
import gwlensing.inference

def main():
    '''
    Main function - takes the user generated ini to compile condor submission and DAG files
    in order to perform the analysis run
    '''

    #Instantiate the Configuration Parser
    config = ConfigParser()

    #If user hasn't given a useable ini raise exception
    if not os.path.isfile(sys.argv[1]):
        raise IOError("Input ini file not given!")

    #Check that the Configuration Parser can read the ini file
    try:
        config.read(sys.argv[1])
    except IOError:
        print("Ini file unreadable!")

    #Create the Outdir Directory and the Data and Submit subfolders
    outdir = config.get("bilby_setup", "outdir")
    data_subfolder = outdir + "/data"
    submit_subfolder = outdir + "/submit"

    for folder in (outdir, data_subfolder, submit_subfolder):
        if not os.path.isdir(folder):
            os.mkdir(folder)

    #Get Injection Parameters and Convert to Floats
    injection_parameters = config._sections["base_waveform_injection_parameters"].copy()
    waveform_arguments = config._sections["waveform_arguments"].copy()

    for key in ["reference_frequency", "minimum_frequency"]:
        waveform_arguments[key] = float(waveform_arguments[key])

    injection_parameters.update((key, float(value)) for key, value in injection_parameters.items())

    #Add the Chirp Mass and Mass Ratio values to the injection parameters
    injection_parameters["chirp_mass"] = bilby.gw.conversion.component_masses_to_chirp_mass(
            injection_parameters["mass_1"], injection_parameters["mass_2"])
    injection_parameters["mass_ratio"] = bilby.gw.conversion.component_masses_to_mass_ratio(
            injection_parameters["mass_1"], injection_parameters["mass_2"])

    #Get the Dimensionless Frequency and Impact Parameters
    w_array_file, y_array_file = gwlensing.lensing.utils.wy_handler(config)

    #Get the Amplification Factor Submit File
    amp_fac_real_file, amp_fac_imag_file = gwlensing.lensing.utils.amp_fac_handler(
            config, w_array_file, y_array_file, mode="pipe")

    #Add Files to the Waveform Arguments
    waveform_arguments["w_array_file"] = w_array_file
    waveform_arguments["y_array_file"] = y_array_file
    waveform_arguments["amp_fac_real_file"] = amp_fac_real_file
    waveform_arguments["amp_fac_imag_file"] = amp_fac_imag_file

    #If Injecting, Generate the Injection File
    if config.getboolean("bilby_pipe_settings", "injection"):
        inject_file = gwlensing.lensing.utils.gen_inject_file(config, injection_parameters)
    else:
        inject_file = None

    #If user specifies, perform unlensed preparation run
    if config.getboolean("data_settings", "create_unlensed_prep_run"):
        #Generate unlensed run ini
        unlensed_ini = gwlensing.lensing.utils.gen_bilby_pipe_ini(
                config=config, inject_file=inject_file,
                waveform_arguments=waveform_arguments,
                mode="unlensed")

        #Perform Initial bilby-pipe run
        subprocess.run(["bilby_pipe", unlensed_ini], check=True)
        subprocess.run(["mv", unlensed_ini, outdir+"/"+unlensed_ini], check=True)

        #TODO: MODIFY PRIORS BASED ON UNLENSED RUN

    #Generate Lensed Run Ini
    lensed_ini = gwlensing.lensing.utils.gen_bilby_pipe_ini(
            config=config, inject_file=inject_file,
            waveform_arguments=waveform_arguments,
            mode="lensed")

    #Perform Initial bilby-pipe run
    subprocess.run(["bilby_pipe", lensed_ini], check=True)
    subprocess.run(["mv", lensed_ini, outdir+"/"+lensed_ini], check=True)

    #Generate Final Overarching DAG
    dag_file = gwlensing.lensing.utils.gen_overarch_dag(config)

    #Message User with Submission
    print("To submit, use \n $ condor_submit_dag " + dag_file)
