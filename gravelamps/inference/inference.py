'''
Gravelamps

Local machine lens data generation and analysis

Written by Mick Wright 2021
'''

import os
import sys
import ast

from configparser import ConfigParser

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

    #Get which methodology is being used
    methodology = config.get("analysis_settings", "methodology")

    #Dependent upon methodology selected, construct waveform arguments dictionary
    if methodology == "interpolate":
        #Get the Dimensionless Frequency and Source Position Files
        dim_freq_file, sour_pos_file = gravelamps.inference.helpers.wy_handler(config)

        #Generate the Amplification Factor Files
        amp_fac_real_file, amp_fac_imag_file = gravelamps.inference.helpers.amp_fac_handler(
            config, dim_freq_file, sour_pos_file, mode="local")

        #Actually construct the dictionary
        analysis_waveform_arguments = gravelamps.inference.helpers.construct_waveform_arguments(
            config, "analysis", dim_freq_file=dim_freq_file, sour_pos_file=sour_pos_file,
            amp_fac_real_file=amp_fac_real_file, amp_fac_imag_file=amp_fac_imag_file)

    elif methodology in ("direct-nonnfw", "direct-nfw"):
        analysis_waveform_arguments = gravelamps.inference.helpers.construct_waveform_arguments(
            config, "analysis")

    #Load in the injection methodolgy to see if it is different
    injection_methodology = config.get("injection_settings", "methodology")

    if injection_methodology is None:
        waveform_arguments = analysis_waveform_arguments.copy()
    else:
        waveform_arguments = gravelamps.inference.helpers.construct_waveform_arguments(
            config, "data")

    #Load in the Priors
    trigger_time = config.getfloat("analysis_settings", "trigger_time")
    duration = config.getfloat("analysis_settings", "duration")
    gps_start_time = trigger_time + 2 - duration
    psd_start_time = trigger_time - duration

    prior_file = config.get("analysis_settings", "prior_file")
    priors = bilby.core.prior.PriorDict(prior_file)

    #Get the Waveform Generator and Frequency Domain Source Model
    lensed_waveform_generator_class = config.get(
        "analysis_settings", "lensed_waveform_generator_class")
    lensed_frequency_domain_source_model = config.get(
        "analysis_settings", "lensed_frequency_domain_source_model")

    lensed_waveform_generator_class, lensed_frequency_domain_source_model = (
        gravelamps.inference.helpers.wfgen_fd_source(
            lensed_waveform_generator_class, lensed_frequency_domain_source_model))

    #Read in Sampling Frequency
    sampling_frequency = config.getfloat("analysis_settings", "sampling_frequency")

    #Set up the Lensed Waveform Generator and Frequency Domain Source Model
    lensed_waveform_generator_class = config.get(
        "analysis_settings", "lensed_waveform_generator_class")
    lensed_frequency_domain_source_model = config.get(
        "analysis_settings", "lensed_frequency_domain_source_model")

    lensed_waveform_generator_class, lensed_frequency_domain_source_model = (
        gravelamps.inference.helpers.wfgen_fd_source(
            lensed_waveform_generator_class, lensed_frequency_domain_source_model))

    #Generate Lensed Waveform
    lensed_waveform_generator = lensed_waveform_generator_class(
        duration=duration, sampling_frequency=sampling_frequency,
        frequency_domain_source_model=lensed_frequency_domain_source_model,
        parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
        waveform_arguments=waveform_arguments)

    #Generate Analysis Waveform
    analysis_waveform_generator = lensed_waveform_generator_class(
        duration=duration, sampling_frequency=sampling_frequency,
        frequency_domain_source_model=lensed_frequency_domain_source_model,
        parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
        waveform_arguments=analysis_waveform_arguments)

    #Set up the Interferometers
    interferometer_list = config.get(
        "analysis_settings", "interferometers").replace(" ","").split(",")

    #Generate the Interferometer Data
    if config.getboolean("injection_settings", "injection"):
        #Create the proper list
        interferometers = bilby.gw.detector.InterferometerList(interferometer_list)

        #Read in the Injection Parameters
        injection_parameters = config._sections["injection_parameters"].copy()
        injection_parameters.update(
            (key,float(value)) for key, value in injection_parameters.items())

        #Add the Chirp Mass and Mass Ratio values to the injection parameters
        injection_parameters["chirp_mass"] = bilby.gw.conversion.component_masses_to_chirp_mass(
            injection_parameters["mass_1"], injection_parameters["mass_2"])
        injection_parameters["mass_ratio"] = bilby.gw.conversion.component_masses_to_mass_ratio(
            injection_parameters["mass_1"], injection_parameters["mass_2"])
        injection_parameters["geocent_time"] = trigger_time

        #Inject Signal into Interferometer
        interferometers.set_strain_data_from_power_spectral_densities(
            sampling_frequency=sampling_frequency, duration=duration, start_time=gps_start_time)
        interferometers.inject_signal(waveform_generator=lensed_waveform_generator,
                                      parameters=injection_parameters)

    else:
        #Create Empty Interferometer List
        interferometers = bilby.gw.detector.InterferometerList([])

        #Load in the Channels
        channels = config.get("event_settings","channels").replace(" ","").split(",")

        #Load in the PSD Files
        psd_files = config.get("event_settings", "psd_files")
        psd_dict = ast.literal_eval(psd_files)
        psd_duration = config.getfloat("event_settings", "psd_duration")

        #Load in the Calibration Files
        calibration_files = config.get("event_settings", "calibration_files")
        calibration_dict = ast.literal_eval(calibration_files)

        #For each channel, load in the event data and set up the calibration
        for channel in channels:
            interferometer = bilby.gw.detector.load_data_by_channel_name(
                channel_name=channel, start_time=gps_start_time, psd_start_time=psd_start_time,
                segment_duration=duration, psd_duration=psd_duration,
                sampling_frequency=sampling_frequency)
            interferometer.minimum_frequency = analysis_waveform_arguments["minimum_frequency"]
            interferometer.maximum_frequency = analysis_waveform_arguments["maximum_frequency"]
            interferometer.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(
                psd_file=psd_dict[interferometer.name])
            interferometer.calibration_model = bilby.gw.calibration.CubicSpline(
                prefix="recalib_" + interferometer.name + "_",
                minimum_frequency=analysis_waveform_arguments["minimum_frequency"],
                maximum_frequency=analysis_waveform_arguments["maximum_frequency"], n_points=10)

            calibration_prior = bilby.gw.prior.CalibrationPriorDict().from_envelope_file(
                calibration_dict[interferometer.name],
                minimum_frequency=analysis_waveform_arguments["minimum_frequency"],
                maximum_frequency=analysis_waveform_arguments["maximum_frequency"],
                n_nodes=10, label=interferometer.name)
            priors.update(calibration_prior)

    #Plot the Interferometer Data
    interferometers.plot_data(outdir=data_subdirectory, label=label)

    #Handle the Sampler Settings
    sampler = config.get("analysis_settings", "sampler")
    sampler_kwargs_dict = config._sections["sampler_kwargs"].copy()

    for key, value in sampler_kwargs_dict.items():
        sampler_kwargs_dict[key] = int(value)

    #If wanting an unlensed analysis run, do this:
    if config.getboolean("unlensed_analysis_settings", "unlensed_analysis_run"):
        #Get the unlensed waveform generator and frequency domain source model
        unlensed_waveform_generator_class = config.get(
            "unlensed_analysis_settings", "unlensed_waveform_generator_class")
        unlensed_frequency_domain_source_model = config.get(
            "unlensed_analysis_settings", "unlensed_frequency_domain_source_model")

        unlensed_waveform_generator_class, unlensed_frequency_domain_source_model = (
            gravelamps.inference.helpers.wfgen_fd_source(
                unlensed_waveform_generator_class, unlensed_frequency_domain_source_model))

        #Generate the Unlensed Waveform
        unlensed_waveform_generator = unlensed_waveform_generator_class(
            duration=duration, sampling_frequency=sampling_frequency,
            frequency_domain_source_model=unlensed_frequency_domain_source_model,
            parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
            waveform_arguments=waveform_arguments)

        #Generate the Unlensed Likelihoods
        unlensed_likelihood = bilby.gw.GravitationalWaveTransient(
            interferometers=interferometers, waveform_generator=unlensed_waveform_generator)

        #Fix the Lens Priors for Unlensed Case
        unlensed_priors = priors

        for parameter in ["lens_mass", "source_position", "lens_fractional_distance"]:
            unlensed_priors[parameter] = 0

        #Perform the Unlensed Run
        if config.getboolean("injection_settings", "injection"):
            result_unlensed = bilby.run_sampler(
                likelihood=unlensed_likelihood, priors=unlensed_priors, outdir=outdir,
                label=label+"_unlensed", sampler=sampler,
                injection_parameters=injection_parameters, **sampler_kwargs_dict)
        else:
            result_unlensed = bilby.run_sampler(
                likelihood=unlensed_likelihood, priors=unlensed_priors, outdir=outdir,
                label=label+"_unlensed", sampler=sampler, **sampler_kwargs_dict)

        if config.getboolean("analysis_settings", "plot_corner"):
            result_unlensed.plot_corner()

        #TODO: Handle resultant posteriors

    #Generate Lensed Likelihood
    lensed_likelihood = bilby.gw.GravitationalWaveTransient(
        interferometers=interferometers, waveform_generator=analysis_waveform_generator)

    #Perform the Lensed Run
    if config.getboolean("injection_settings", "injection"):
        result_lensed = bilby.run_sampler(
            likelihood=lensed_likelihood, priors=priors, outdir=outdir, label=label,
            injection_parameters=injection_parameters, sampler=sampler, **sampler_kwargs_dict)
    else:
        result_lensed = bilby.run_sampler(
            likelihood=lensed_likelihood, priors=priors, outdir=outdir, label=label,
            sampler=sampler, **sampler_kwargs_dict)

    if config.getboolean("analysis_settings", "plot_corner"):
        result_lensed.plot_corner()
