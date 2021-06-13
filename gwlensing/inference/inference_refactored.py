'''
Gravelamps

Local machine lens data generation and analysis variant
for Simulated Data

Written by Mick Wright 2021
'''

import os
import sys
import importlib

from configparser import ConfigParser

import numpy as np
import bilby

import gwlensing.lensing

def main():
    '''
    Main function - takes the user generated ini and generates the lens data
    and performs an analysis run
    '''

    #Instantiate the Configuration Parser
    config = ConfigParser()

    #If the user hasn't given a useable ini file, raise excpetion
    if not os.path.isfile(sys.argv[1]):
        raise IOError("Input ini file not given!")

    #Check that the Configuration Parser can read the ini file
    try:
        config.read(sys.argv[1])
    except IOError:
        print("Input ini file is unreadable!")

    #Bilby Logging Set Up
    label = config.get("bilby_setup", "label")
    outdir = config.get("bilby_setup", "outdir")

    bilby.core.utils.setup_logger(label=label, outdir=outdir)

    #If the data subdirectory doesn't exist, create it
    data_subdir = config.get("data_settings", "data_subdir")
    data_fullpath = outdir + "/" + data_subdir

    if not os.path.isdir(data_fullpath):
        os.mkdir(data_fullpath)

    #Read in User Parameters for the Bilby Analysis
    duration = config.getfloat("bilby_setup", "duration")
    sampling_frequency = config.getfloat("bilby_setup", "sampling_frequency")

    #Read in the Waveform parameters and Convert to Floats
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

    #Get the Dimensionless Frequency and Impact Parameter Files
    w_array_file, y_array_file = gwlensing.lensing.utils.wy_handler(config)

    #Generate Lensed Interpolator
    #Generate Lensed Waveform
    #Generate Injected Interformeters
    #If unlensed prep run, do this
    #Lensed run
