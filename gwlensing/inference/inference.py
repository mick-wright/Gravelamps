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
    amp_fac_real_file, amp_fac_imag_file = gwlensing.lensing.utils.amp_fac_handler(
            config, w_array_file, y_array_file, mode="local")

    #Generate Lensed Waveform
    lensed_waveform_generator = gwlensing.lensing.Lensed_Waveform_Generator(
            duration=duration, sampling_frequency=sampling_frequency,
            frequency_domian_source_model = gwlensing.lensing.BBH_lensed_waveform,
            parameter_conversion = bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
            waveform_arguments = waveform_arguments)

    #Generate Injected Interformeters
    interferometer_list = config.get("bilby_setup", "detectors").replace(" ","").split(",")
    interferometers = bilby.gw.detector.InterferometerList(interferometer_list)
    interferometers.set_strain_data_from_power_spectral_densities(
            sampling_frequency = sampling_frequency, duration=duration,
            start_time = injection_parameters["geocent_time"] - 3)
    interferometers.inject_signal(waveform_generator=lensed_waveform_generator,
            parameters = injection_parameters)

    #Loading in the Prior File, Fixing Any Parameters Specified
    prior_file = config.get("prior_settings", "prior_file")
    prior_fix_list = config.get("prior_settings", "parameters_to_fix").replace(" ","").split(",")

    priors = bilby.core.prior.PriorDict(prior_file)
    for parameter in prior_fix_list:
        priors[parameter] = injection_parameters[parameter]

    #Sampler Settings
    sampler = config.get("bilby_setup", "sampler")
    plot_corner = config.getboolean("bilby_setup", "plot_corner")
    sampler_kwargs_dict = config._sections["sampler_kwargs"].copy()

    #If unlensed prep run, do this
    if config.getboolean("data_settings", "create_unlensed_prep_run"):
        #Generate the Unlensed Waveform
        unlensed_waveform_generator = bilby.gw.WaveformGenerator(
                duration = duration, sampling_frequency = sampling_frequency,
                frequency_domain_source_model = bilby.gw.source.lal_binary_black_hole,
                parameter_conversion = bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
                waveform_arguments = waveform_arguments)

        #Get the Unlensed Likelihoods
        unlensed_likelihood = bilby.gw.GravitationalWaveTransient(
                interferometers = interferometers,
                waveform_generator = unlensed_waveform_generator)

        #Fix Lens Priors for Unlensed Case
        unlensed_priors = priors

        for parameter in ["lens_mass", "impact_parameter", "lens_fractional_distance"]:
            unlensed_priors[parameter] = 0

        #Perform the Unlensed Run
        result_unlensed = bilby.run_sampler(
                likelihood = unlensed_likelihood, priors = unlensed_priors,
                outdir = outdir, label = label+"_unlensed", sampler = sampler,
                plot = plot_corner, injection_parameters = injection_parameters,
                **sampler_kwargs_dict)

        #TODO: Get either the maximum posterior parameter values, or the full posteriors to use as priors

    #Lensed run
    #Get the Likelihood
    lensed_likelihood = bilby.gw.GravitationalWaveTransient(
            interferometers = interferometers,
            waveform_generator = lensed_waveform_generator)

    #Run the Sampler
    result_lensed = bilby.run_sampler(
            likelihood = lensed_likelihood, priors = priors,
            outdir = outdir, label = label, sampler = sampler,
            plot = plot_corner, injection_parameters = injection_parameters,
            **sampler_kwargs_dict)


