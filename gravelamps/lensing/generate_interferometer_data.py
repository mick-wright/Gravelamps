'''
Gravelamps

Program generates interferometer data - to simulate real event handling.
This allows simulating one type of lens given another.

N.B: REQUIRES THAT THE LENSING DATA BE GIVEN IN THE INJECTION SETTINGS

Written by Mick Wright 2021
'''

import os
import sys

from configparser import ConfigParser

import bilby

import gravelamps.lensing
import gravelamps.inference

def main():
    '''
    Main Function - takes in the user generated INI and generates the interferometer data
    based upon the injection settings. Requires that the lensing data already be generated
    and given in the settings - this makes it cluster friendly as well as locally runable.
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
        print("Input INI file cannot be read!")

    #Bilby Logging Set Up
    label = config.get("output_settings", "label")
    outdir = config.get("output_settings", "outdir")

    #Get Data Subdirectory
    data_subdirectory = outdir + "/data"

    if not os.path.isdir(data_subdirectory):
        os.mkdir(data_subdirectory)

    bilby.core.utils.setup_logger(label=label, outdir=outdir)

    #Get the Waveform Generator and Frequency Domain Source Model
    lensed_waveform_generator_class = config.get(
        "analysis_settings", "lensed_waveform_generator_class")
    lensed_frequency_domain_source_model = config.get(
        "analysis_settings", "lensed_frequency_domain_source_model")

    lensed_waveform_generator_class, lensed_frequency_domain_source_model = (
        gravelamps.inference.helpers.wfgen_fd_source(
            lensed_waveform_generator_class, lensed_frequency_domain_source_model))

    #Load in duration, and start times
    trigger_time = config.getfloat("analysis_settings", "trigger_time")
    duration = config.getfloat("analysis_settings", "duration")
    gps_start_time = trigger_time + 2 - duration

    #Construct the Waveform Arguments Dictionary
    dim_freq_file = os.path.abspath(
        config.get("injection_settings", "dimensionless_frequency_file"))
    sour_pos_file = os.path.abspath(
        config.get("injection_settings", "source_position_file"))
    amp_fac_real_file = os.path.abspath(
        config.get("injection_settings", "amplification_factor_real_file"))
    amp_fac_imag_file = os.path.abspath(
        config.get("injection_settings", "amplification_factor_imag_file"))

    waveform_approximant = config.get("analysis_settings", "waveform_approximant")
    minimum_frequency = config.get("analysis_settings", "minimum_frequency")
    maximum_frequency = config.get("analysis_settings", "maximum_frequency")
    reference_frequency = config.get("analysis_settings", "reference_frequency")

    waveform_arguments = {}

    waveform_arguments["waveform_approximant"] = waveform_approximant
    waveform_arguments["minimum_frequency"] = minimum_frequency
    waveform_arguments["maximum_frequency"] = maximum_frequency
    waveform_arguments["reference_frequency"] = reference_frequency
    waveform_arguments["dim_freq_file"] = dim_freq_file
    waveform_arguments["sour_pos_file"] = sour_pos_file
    waveform_arguments["amp_fac_real_file"] = amp_fac_real_file
    waveform_arguments["amp_fac_imag_file"] = amp_fac_imag_file

    #Read in the Sampling Frequnecy
    sampling_frequency = config.getfloat("analysis_settings", "sampling_frequency")

    #Create Waveform Generator
    lensed_waveform_generator = lensed_waveform_generator_class(
        duration=duration, sampling_frequency=sampling_frequency,
        frequency_domain_source_model=lensed_frequency_domain_source_model,
        parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
        waveform_arguments=waveform_arguments)

    #Set up Interferomters
    interferometer_list = config.get(
        "analysis_settings", "interferometers").replace(" ", "").split(",")

    #Generate the Interferometer Data
    interferometers = bilby.gw.detector.InterferometerList(interferometer_list)

    injection_parameters = config._sections["injection_parameters"].copy()
    injection_parameters.update(
        (key, float(value)) for key, value in injection_parameters.items())

    #Add the Chirp Mass and Mass Ratio Values to the Injection Parameters
    injection_parameters["chirp_mass"] = bilby.gw.conversion.component_masses_to_chirp_mass(
        injection_parameters["mass_1"], injection_parameters["mass_2"])
    injection_parameters["mass_ratio"] = bilby.gw.conversion.component_masses_to_mass_ratio(
        injection_parameters["mass_1"], injection_parameters["mass_2"])
    injection_parameters["geocent_time"] = trigger_time

    #Inject Signal Into Interferometer
    interferometers.set_strain_data_from_power_spectral_densities(
        sampling_frequency=sampling_frequency, duration=duration, start_time=gps_start_time)
    interferometers.inject_signal(waveform_generator=lensed_waveform_generator,
                                  parameters=injection_parameters)

    #Save Data from Interferometers to data subdirectory
    for ifo in interferometers:
        ifo.savedata(data_subdirectory)
