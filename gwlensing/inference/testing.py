import sys
import os
import numpy as np
import bilby
import gwlensing.lensing

from configparser import ConfigParser

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

    #Bilby Setup
    duration = config.getfloat("bilby_setup","duration")
    sampling_frequency = config.getfloat("bilby_setup","sampling_frequency")

    #Bilby Logging Set Up
    bilby.core.utils.setup_logger(label=config.get("bilby_setup","label"), outdir=config.get("bilby_setup","outdir"))

    #Injection Parameters for Base Waveform
    injection_parameters = config._sections["base_waveform_injection_parameters"]
    waveform_arguments = config._sections["waveform_arguments"]

    #Convert to floats
    waveform_arguments["reference_frequency"] = float(waveform_arguments["reference_frequency"])
    waveform_arguments["minimum_frequency"] = float(waveform_arguments["minimum_frequency"])
    injection_parameters.update((k,float(v)) for k, v in injection_parameters.items())

    #Generate Lens Interpolator
    if config.get("optional_input","w_array_file") == "None":
        #Calculate Redshifted Lens Mass
        lens_distance = injection_parameters["luminosity_distance"]*injection_parameters["lens_fractional_distance"]

        lens_redshift = bilby.gw.conversion.luminosity_distance_to_redshift(lens_distance)

        redshifted_lens_mass = gwlensing.lensing.utils.natural_mass(injection_parameters["lens_mass"]*(1+lens_redshift))

        w_array = gwlensing.lensing.utils.generate_dimensionless_frequency_array(1250,redshifted_lens_mass)
    else:
        try:
            w_array = np.loadtxt(config.get("optional_input","w_array_file"))
        except IOError:
            print("W Array File Couldn't Be Read!")

    if config.get("optional_input","y_array_file") == "None":
        y_array = np.linspace(0.01, 1.10, 20)
    else:
        try:
            y_array = np.loadtxt(config.get("optional_input","y_array_file"))
        except IOError:
            print("Y Array file Couldn't Be Read")

    if config.get("optional_input","amp_fac_file") == "None":
        f_matrix = gwlensing.lensing.point_lens.generate_amplification_factor_matrix(w_array, y_array)
    else:
        try:
            f_matrix = np.loadtxt(config.get("optional_input","amp_fac_file"), dtype=complex)
        except IOError:
            print("Amplification Factor file couldn't be read!")

    lens_interpolator = gwlensing.lensing.utils.generate_interpolator(w_array, y_array, f_matrix)

    waveform_arguments["interpolator"] = lens_interpolator

    wfgen = bilby.gw.WaveformGenerator(duration=duration, sampling_frequency=sampling_frequency, frequency_domain_source_model=gwlensing.lensing.BBH_lensed_waveform, parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters, waveform_arguments=waveform_arguments)

    interferometer_list = config.get("bilby_setup","detectors").replace(" ","").split(",")

    ifos = bilby.gw.detector.InterferometerList(interferometer_list)
    ifos.set_strain_data_from_power_spectral_densities(sampling_frequency=sampling_frequency,duration=duration, start_time=injection_parameters["geocent_time"]-3)
    ifos.inject_signal(waveform_generator=wfgen,parameters=injection_parameters)

    likelihood = bilby.gw.GravitationalWaveTransient(interferometers=ifos, waveform_generator=wfgen)

    priors = bilby.core.prior.PriorDict(config.get("prior_settings","prior_file"))

    prior_fix_list = config.get("prior_settings","parameters_to_fix").replace(" ","").split(",")

    for parameter in prior_fix_list:
        priors[parameter] = injection_parameters[parameter]

    result = bilby.run_sampler(likelihood=likelihood, priors=priors, sampler=config.get("sampler_settings","sampler"), nlive=config.getint("sampler_settings","nlive"), npool=config.getint("sampler_settings","npool"),nact=config.getint("sampler_settings","nact"),nparallel=config.getint("sampler_settings","nparallel"),injection_parameters=injection_parameters, outdir=config.get("bilby_setup","outdir"), label=config.get("bilby_setup","label"))

    if config.getbool("bilby_setup","make_corner") == True:
        result.plot_corner()
