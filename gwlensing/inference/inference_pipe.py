import sys
import os
import numpy as np
import bilby_pipe
import bilby
import gwlensing.lensing
import gwlensing.inference

from configparser import ConfigParser

def main():
    #Instantiate the Configuration Parser
    config = ConfigParser() 

    #If user hasn't given a useable ini raise exception
    if not os.path.isfile(sys.argv[1]):
        raise IOError("Input ini file not given!")

    #Check that the Configuration Parser can read hte ini file
    try:
        config.read(sys.argv[1])
    except IOError:
        print("Ini file unreadable!") 

    #Get in Injection Parameters and Convert to Floats 
    injection_parameters = config._sections["base_waveform_injection_parameters"].copy()
    waveform_arguments = config._sections["waveform_arguments"].copy()

    waveform_arguments["reference_frequency"] = float(waveform_arguments["reference_frequency"])
    waveform_arguments["minimum_frequency"] = float(waveform_arguments["minimum_frequency"]) 
    injection_parameters.update((key, float(value)) for key, value in injection_parameters.items())

    injection_parameters["chirp_mass"] = bilby.gw.conversion.component_masses_to_chirp_mass(injection_parameters["mass_1"], injection_parameters["mass_2"])

    #Get/Generate Dimensionless Frequency and Impact Parameter Files
    w_array_file, y_array_file = gwlensing.lensing.utils.wyhandler(config, injection_parameters)

    #Get Lens Model
    lens_model = config.get("lens_settings", "lens_model")
    #Get Additional Lensing Parameters if neeeded
    additional_lens_parameters = gwlensing.lensing.utils.get_additional_parameters(config) 

    #Create the Generate Amplification Factor submit file - if needed
    amp_fac_real_file, amp_fac_imag_file = gwlensing.lensing.utils.ampfachandler(config, injection_parameters, w_array_file, y_array_file, lens_model, additional_lens_parameters=[], mode="pipe") 
