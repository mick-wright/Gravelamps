'''
Gravelamps

Local lens data generation only

Written by Mick Wright 2021
'''

import os
import sys

from configparser import ConfigParser

import gravelamps.lensing
import gravelamps.inference

def main():
    '''
    Main function - takes the user generated ini and generates the lens data
    only
    '''

    #Instantiate the Configuration Parser
    config = ConfigParser()

    #If the user hasn't given a useable INI file, raise exception
    if not os.path.isfile(sys.argv[1]):
        raise IOError("Input INI file not given!")

    #Check that the Configuration Parser can read the INI
    try:
        config.read(sys.argv[1])
    except IOError:
        print("Input INI cannot be read!")

    #Generate the Lens Data
    output_directory = config.get("output_settings", "outdir")
    data_subdirectory = output_directory + "/data"

    #If the Data subdirectory doesn't exist make it, along with output directory
    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)
    if not os.path.isdir(data_subdirectory):
        os.mkdir(data_subdirectory)

    #Get the Dimensionless Frequency and Source Position Files
    dim_freq_file, sour_pos_file = gravelamps.inference.helpers.wy_handler(config)

    #Generate the Amplification Factor Files
    amp_fac_real_file, amp_fac_imag_file = gravelamps.inference.helpers.amp_fac_handler(
        config, dim_freq_file, sour_pos_file, mode="local")
