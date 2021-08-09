'''
Gravelamps

Local machine lens data generation and analysis

Written by Mick Wright 2021
'''

import os
import sys

from configparser import ConfigParser

import numpy as np
import bilby

import gravelamps.lensing
import gravelamps.inference

def main():
    '''
    Main function - takes the user generated ini and generates the lens data
    and performs an analysis run
    '''

    #Instantiate the Configuration Parser
    config = ConfigParser()

    #If the user hasn't given a useable ini file, raise exception
    if not os.path.isfile(sys.argv[1]):
        raise IOError("Input INI file not given!")

    #Check that the Configuration Parser can read the INI
    try:
        config.read(sys.argv[1])
    except IOError:
        print("Input ini file cannot be read!")

    #Bilby Logging set up
    label = config.get("output_settings", "label")
    outdir = config.get("output_settings", "outdir")

    bilby.core.utils.setup_logger(label=label, outdir=outdir)

    #Generate the Lensing Interpolator
    #If the data subdirectory doesn't already exist, create it
    data_subdirectory = outdir + "/data"

    if not os.path.isdir(data_subdirectory):
        os.mkdir(data_subdirectory)

    #Get the Dimensionless Frequency and Source Position Files
    dim_freq_file, sour_pos_file = gravelamps.inference.helpers.wy_handler(config)

    #Generate the Amplification Factor Files
    amp_fac_real_file, amp_fac_imag_file = gravelamps.inference.helpers.amp_fac_handler(
        config, dim_freq_file, sour_pos_file, mode="local")
