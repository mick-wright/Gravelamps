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
		print("Ini file unreadable!") 

	#Bilby Logging Set Up 
	label = config.get("bilby_setup","label")
	outdir = config.get("bilby_setup","outdir")

	bilby.core.utils.setup_logger(label=label, outdir=outdir)
	
	data_subdir = config.get("data_settings","data_subdir")
	if not os.path.isdir(outdir+"/"+data_subdir):
		os.mkdir(outdir+"/"+data_subdir) 

	#Read in User Parameters for Bilby Analysis
	duration = config.getfloat("bilby_setup","duration")
	sampling_frequency = config.getfloat("bilby_setup","sampling_frequency") 

	#Read in Waveform Parameters and convert to floats
	injection_parameters = config._sections["base_waveform_injection_parameters"].copy()
	waveform_arguments = config._sections["waveform_arguments"].copy()

	waveform_arguments["reference_frequency"] = float(waveform_arguments["reference_frequency"]) 
	waveform_arguments["minimum_frequency"] = float(waveform_arguments["minimum_frequency"]) 
	injection_parameters.update((key, float(value)) for key, value in injection_parameters.items()) 

	#Get the Dimensionless Frequency and Impact Parameter Files 
	w_array_file, y_array_file = gwlensing.lensing.utils.wyhandler(config, injection_parameters)

	#Get the Amplification Factor Files
	amp_fac_real_file, amp_fac_imag_file = gwlensing.lensing.utils.ampfachandler(config,injection_parameters, w_array_file, y_array_file) 

	#Generate the Lens Interpolator Function and add it to the waveform arguments 
	lens_interpolator = gwlensing.lensing.utils.generate_interpolator(w_array_file, y_array_file, amp_fac_real_file, amp_fac_imag_file)
	waveform_arguments["interpolator"] = lens_interpolator 

	#Generate the Lensed Waveform
	wfgen = bilby.gw.WaveformGenerator(duration=duration, sampling_frequency=sampling_frequency, frequency_domain_source_model=gwlensing.lensing.BBH_lensed_waveform, parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters, waveform_arguments=waveform_arguments)

	#Create Interferometers and inject lensed signal
	interferometer_list = config.get("bilby_setup","detectors").replace(" ","").split(",")
	ifos = bilby.gw.detector.InterferometerList(interferometer_list) 
	ifos.set_strain_data_from_power_spectral_densities(sampling_frequency=sampling_frequency,duration=duration, start_time=injection_parameters["geocent_time"]-3)
	ifos.inject_signal(waveform_generator=wfgen, parameters=injection_parameters) 

	#Likelihood 
	likelihood = bilby.gw.GravitationalWaveTransient(interferometers=ifos, waveform_generator=wfgen)
	
	#Load in Prior File and Fix Specified Parameters 
	priors = bilby.core.prior.PriorDict(config.get("prior_settings","prior_file")) 
	prior_fix_list = config.get("prior_settings","parameters_to_fix").replace(" ","").split(",") 
	for parameter in prior_fix_list:
		priors[parameter] = injection_parameters[parameter] 

	#Run Sampler, if user has specified generate a corner plot
	result = bilby.run_sampler(likelihood=likelihood, priors=priors, sampler=config.get("sampler_settings","sampler"), nlive=config.getint("sampler_settings","nlive"), npool=config.getint("sampler_settings","npool"),nact=config.getint("sampler_settings","nact"),nparallel=config.getint("sampler_settings","nparallel"),injection_parameters=injection_parameters, outdir=config.get("bilby_setup","outdir"), label=config.get("bilby_setup","label"))
	if config.getboolean("bilby_setup","make_corner") == True:
		result.plot_croner() 

	#Generate a "Treat Unlensed" if user has specified they want one
	if config.getboolean("treat_unlensed_settings","create_treat_unlensed") == True: 
		wfgen_unlensed = bilby.gw.WaveformGenerator(duration=duration, sampling_frequency=sampling_frequency, frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole, parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters, waveform_arguments=waveform_arguments)

		likelihood_unlensed = bilby.gw.GravitationalWaveTransient(interferometers=ifos, waveform_generator=wfgen_unlensed)

		#For Unlensed Priors, fix all lens parameters 
		priors_unlensed = priors
		for parameter in ["lens_mass","impact_parameter","lens_fractional_distance"]:
			priors_unlensed[parameter] = injection_parameters[parameter] 

		result_unlsned = bilby.run_sampler(likelihood=likelihood_unlensed, priors=priors_unlensed, sampler=config.get("sampler_settings","sampler"), nlive=config.getint("sampler_settings","nlive"), npool=config.getint("sampler_settings","npool"),nact=config.getint("sampler_settings","nact"),nparallel=config.getint("sampler_settings","nparallel"),injection_parameters=injection_parameters, outdir=config.get("treat_unlensed_settings","unlensed_outdir"), label=config.get("treat_unlensed_settings","unlensed_label"))

		#Generate a "Treat Unlensed" Corner Plot if user desires
		if config.getboolean("treat_unlensed_settings","unlensed_make_corner") == True:
			result_unlensed.plot_corner() 
