'''
Gravelamps

Clusterised lens data generation only
Assumes cluster is using HTCondor Scheduler

Written by Mick Wright 2021
'''

import os
import sys

from configparser import ConfigParser

import gravelamps.lensing
import gravelamps.inference

def main():
    '''
    Main function - takes the user generated INI and writes the lens generation
    submission file only
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
        print("Input INI cannot be read!")

    #Create the Output Directory and Data and Submit Subdirectories
    outdir = config.get("output_settings", "outdir")
    data_subdirectory = outdir + "/data"
    submit_subdirectory = outdir + "/submit"

    for folder in (outdir, data_subdirectory, submit_subdirectory):
        if not os.path.isdir(folder):
            os.mkdir(folder)

    #Generate the Lensing Data
    #Get the Dimensionless Frequency and Source Position Files
    dim_freq_file, sour_pos_file = gravelamps.inference.helpers.wy_handler(config)

    #Generate the Amplification Factor Submit File
    amp_fac_real_file, amp_fac_imag_file = gravelamps.inference.helpers.amp_fac_handler(
        config, dim_freq_file, sour_pos_file, mode="pipe")

    #Give user submission message
    print("Submit file located in the submit subdirectory: {0}".format(submit_subdirectory))
    print("Results upon completion will be found in the data subdirectory: {0}".format(data_subdirectory)) 

