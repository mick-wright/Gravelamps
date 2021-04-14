import sys
import os
import numpy as np
import bilby_pipe
import bilby
import gwlensing.lensing
import gwlensing.inference
import subprocess

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

    #Create the Outdir Directory and the Data and Submit subfolders
    outdir = config.get("bilby_setup", "outdir")

    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    if not os.path.isdir(outdir+"/data"):
        os.mkdir(outdir+"/data")
    if not os.path.isdir(outdir+"/submit"):
        os.mkdir(outdir+"/submit") 

    #Get in Injection Parameters and Convert to Floats 
    injection_parameters = config._sections["base_waveform_injection_parameters"].copy()
    waveform_arguments = config._sections["waveform_arguments"].copy()

    waveform_arguments["reference_frequency"] = float(waveform_arguments["reference_frequency"])
    waveform_arguments["minimum_frequency"] = float(waveform_arguments["minimum_frequency"]) 
    injection_parameters.update((key, float(value)) for key, value in injection_parameters.items())

    #Get/Generate Dimensionless Frequency and Impact Parameter Files
    w_array_file, y_array_file = gwlensing.lensing.utils.wyhandler(config, injection_parameters)

    waveform_arguments["w_array_file"] = os.path.abspath(w_array_file)
    waveform_arguments["y_array_file"] = os.path.abspath(y_array_file)

    #Get Lens Model
    lens_model = config.get("lens_settings", "lens_model")
    #Get Additional Lensing Parameters if neeeded
    additional_lens_parameters = gwlensing.lensing.utils.get_additional_parameters(config) 

    #Create the Generate Amplification Factor submit file - if needed
    amp_fac_real_file, amp_fac_imag_file = gwlensing.lensing.utils.ampfachandler(config, injection_parameters, w_array_file, y_array_file, lens_model, additional_lens_parameters, mode="pipe")

    waveform_arguments["amp_fac_real_file"] = os.path.abspath(amp_fac_real_file)
    waveform_arguments["amp_fac_imag_file"] = os.path.abspath(amp_fac_imag_file)

    #Generate Injection File
    inject_file = gwlensing.lensing.utils.gen_inject_file(config, outdir, injection_parameters)

    #Create the bilby-pipe ini file
    gwlensing.lensing.utils.gen_bilby_pipe_ini(config, outdir, inject_file, waveform_arguments) 

    #Run Initial bilby-pipe run 
    subprocess.run(["bilby_pipe", config.get("bilby_setup","label")+"_bilby_pipe.ini"])
    subprocess.run(["rm", config.get("bilby_setup","label")+"_bilby_pipe.ini"]) 

    #Generate final overarching dag
    gwlensing.lensing.utils.gen_overarch_dag(config, outdir)

    #Message User with Submission
    print("To submit, use \n $ condor_submit_dag " + outdir+"/submit/dag_"+config.get("bilby_setup","label")+"_overarch.submit")
