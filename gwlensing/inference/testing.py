import sys
import os
import bilby

from configparser import ConfigParser

from ..lensing import point_lens as pl

def main():
    #Instantiate the Configuration Parser
    config = ConfigParser()

    #If user hasn't given a useable ini file, raise exception
    if not os.path.isfile(sys.argv[1]):
        raise IOError("Input ini file not given")

    #Check that the Configuration Parser can read the ini file
    try:
        config.read(sys.argv[1])
    except IOError:
        print("Input file cannot be read")

    #Bilby Logging Set Up
    bilby.core.utils.setup_logger(label=config.get("bilby_setup","label"), outdir=config.get("bilby_setup","outdir"))

    #Injection Parameters for Base Waveform
    injection_parameters = config._section["base_waveform_injection_parameters"]
    waveform_arguments = config._section["waveform_arguments"]

    #Lens Parameters
    true_impact_parameter = config.getfloat("lens_parameters","impact_parameter")
    true_lens_mass = config.getfloat("lens_parameters","M_lens")
    true_lens_fractional_distance = config.getfloat("lens_parameters","lens_fractional_distance")
