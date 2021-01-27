import numpy as np
import os
import gwlensing.lensing
import warnings
import scipy.interpolate as scint

from configparser import ConfigParser

def generate_dimensionless_frequency_file(config, injection_parameters):
	outdir = config.get("bilby_setup","outdir") 
	data_subdir = outdir+"/"+config.get("data_settings","data_subdir")

	lens_distance = injection_parameters["luminosity_distance"]*injection_parameters["lens_fractional_distance"]
	lens_redshift = bilby.gw.conversion.luminosity_distance_to_redshift(lens_distance)
	redshifted_lens_mass = gwlensing.lensing.natural_mass(lens_mass*(1+lens_redshift))

	minimum_frequency = config.get("waveform_arguments","minimum_frequency")
	maximal_frequency = config.get("data_settings","maximum_frequency")
	changeover_frequency = config.get("data_settings","changeover_frequency")
	below_npoints = config.get("data_settings","below_npoints") 
	above_npoints = config.get("data_settings","above_npoints") 

	lower_frequency_array = np.linspace(minimum_frequency,changeover_frequency,below_npoints)
	higher_frequency_array = np.linspace(changeover_frequency,maximum_frequency,above_npoints+1) 
	higher_frequency_array = np.delete(higher_frequency_array, 0)

	complete_frequency_array = np.concatenate((lower_frequency_array,higher_frequency_array))

	dim_freq_array = gwlensing.lensing.dimensionless_frequency(complete_frequency_array, redshifted_lens_mass) 

	np.savetxt(data_subdir+"/w.dat", dim_freq_array) 

	return(data_subdir+"/w.dat")

def generate_impact_parameter_file(config,injection_parameters):
	outdir = config.get("bilby_setup","outdir")
	data_subdir = outdir+"/"+config.get("data_settings","data_subdir")
	
	y_npoints = config.getfloat("data_settings","y_npoints")
	min_y = config.getfloat("data_settings","min_y")
	max_y = config.getfloat("data_settings","max_y") 

	y_array = np.linspace(min_y, max_y, y_npoints) 

	np.savetxt(data_subdir+"/y.dat", y_array)

	return(data_subdir+"/y.dat") 

def wyhandler(config, injection_parameters):
	outdir = config.get("bilby_setup","outdir") 
	data_subdir = outdir+"/"+config.get("data_settings","data_subdir") 

	if config.get("optional_input","w_array_file") != "None":
		if config.get("optional_input","y_array_file") == "None":
			warnings.warn("y array not found, interpolator may not be accurate in y!")
		w_array_file = config.get("optional_input","w_array_file")
	elif os.path.isfile(data_subdir+"/w.dat"):
		w_array_file = data_subdir+"/w.dat"
	else:
		w_array_file = gwlensing.lensing.utils.generate_dimensionless_frequency_file(config, injection_parameters)

	if config.get("optional_input","y_array_file") != "None":
		if config.get("optional_input","w_array_file") == "None": 
			warnings.warn("w array not found, interpolator may not be accurate in w!")
		y_array_file = config.get("optional_input","y_array_file")
	elif os.path.isfile(data_subdir+"/y.dat"):
		y_array_file = data_subdir+"/y.dat" 
	else:
		y_array_file = gwlensing.lensing.utils.generate_impact_parameter_file(config, injection_parameters)

	if config.getboolean("data_settings","copy_data_files") == True:
		if w_array_file != data_subdir+"/w.dat":
			subprocess.run(["cp", w_array_file, data_subdir+"/w.dat"])
		elif y_array_file != data_subdir+"/y.dat":
			subprocess.run(["cp", y_array_file, data_subdir+"/y.dat"]) 

	return(w_array_file, y_array_file) 

def ampfachandler(config, injection_parameters, w_array_file, y_array_file):
	outdir = config.get("bilby_setup","outdir")
	data_subdir = outdir+"/"+config.get("data_settings","data_subdir") 

	if config.get("optional_input","amp_fac_complex_file") != "None":
		complex_file = config.get("optional_input","amp_fac_complex_file") 
		complex_array = np.loadtxt(complex_file,dtype=complex) 
		real_array = np.real(complex_array)
		imag_array = np.imag(imag_array)
		
		amp_fac_real_file = data_subdir+"/fReal.dat"
		amp_fac_imag_file = data_subdir+"/fImag.dat" 

		np.savetxt(amp_fac_real_file, real_array)
		np.savetxt(amp_fac_imag_file, imag_array) 
	elif config.get("optional_input","amp_fac_real_file") != "None" and config.get("optional_input","amp_fac_imag_file") != "None":
		amp_fac_real_file = config.get("optional_input","amp_fac_real_file")
		amp_fac_imag_file = config.get("optional_input","amp_fac_imag_file") 
	elif os.path.isfile(data_subdir+"/fReal.dat") and os.path.isfile(data_subdir+"/fImag.dat"):
		amp_fac_real_file = data_subdir+"/fReal.dat"
		amp_fac_imag_file = data_subdir+"/fImag.dat" 
	else:
		subprocess.run(["plCalc", w_array_file, y_array_file])
		subprocess.run(["mv", "fReal.dat", "fImag.dat", data_subdir+"/."]) 

		amp_fac_real_file = data_subdir+"/fReal.dat" 
		amp_fac_imag_file = data_subdir+"/fImag.dat" 

	if config.getboolean("data_settings","copy_data_files"):
		if amp_fac_real_file != data_subdir+"/fReal.dat":
			subprocess.run(["cp", amp_fac_real_file, data_subdir+"/fReal.dat"])
		if amp_fac_imag_file != data_subdir+"/fImag.dat": 
			subprocess.run(["cp", amp_fac_imag_file, data_subdir+"/fImag.dat"]) 

	return(amp_fac_real_file, amp_fac_imag_file) 

def generate_interpolator(w_array_file, y_array_file, amp_fac_real_file, amp_fac_imag_file):
	w_array = np.loadtxt(w_array_file) 
	y_array = np.loadtxt(y_array_file)
	amp_fac_real = np.loadtxt(amp_fac_real_file) 
	amp_fac_imag = np.loadtxt(amp_fac_imag_file) 

	amp_interp_real = scint.RectBivariateSpline(w_array, y_array, amp_fac_real, kx=1, ky=1, s=0)
	amp_interp_imag = scint.RectBivariateSpline(w_array, y_array, amp_fac_imag, kx=1, ky=1, s=0) 

	interpolator_func = lambda w, y: (amp_interp_real(w,y)+1j*amp_interp_imag(w,y)).flatten()

	return(interpolator_func) 
