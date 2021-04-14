import numpy as np
import os
import gwlensing.lensing
import warnings
import scipy.interpolate as scint
import bilby
import subprocess

from configparser import ConfigParser

def generate_dimensionless_frequency_file(config, injection_parameters):
	outdir = config.get("bilby_setup","outdir") 
	data_subdir = outdir+"/"+config.get("data_settings","data_subdir")

	lens_distance = injection_parameters["luminosity_distance"]*injection_parameters["lens_fractional_distance"]
	lens_redshift = bilby.gw.conversion.luminosity_distance_to_redshift(lens_distance)
	redshifted_lens_mass = gwlensing.lensing.natural_mass(injection_parameters["lens_mass"]*(1+lens_redshift))

	maximum_frequency = config.getint("data_settings","maximum_frequency")
	changeover_frequency = config.getint("data_settings","changeover_frequency")

	below_npoints = config.getint("data_settings","below_npoints") 
	above_npoints = config.getint("data_settings","above_npoints") 

	lower_frequency_array = np.linspace(0,changeover_frequency,below_npoints)
	higher_frequency_array = np.linspace(changeover_frequency,maximum_frequency,above_npoints+1) 
	higher_frequency_array = np.delete(higher_frequency_array, 0)

	complete_frequency_array = np.concatenate((lower_frequency_array,higher_frequency_array))

	dim_freq_array = gwlensing.lensing.dimensionless_frequency(complete_frequency_array, redshifted_lens_mass) 

	np.savetxt(data_subdir+"/w.dat", dim_freq_array) 

	return(data_subdir+"/w.dat")

def generate_impact_parameter_file(config,injection_parameters):
	outdir = config.get("bilby_setup","outdir")
	data_subdir = outdir+"/"+config.get("data_settings","data_subdir")
	
	y_npoints = config.getint("data_settings","y_npoints")
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

def get_additional_parameters(config):
	lens_model = config.get("lens_settings", "lens_model")

	if lens_model == "pointlens":
		return([])
	elif lens_model == "sislens":
		return([])
	elif lens_model == "nfwlens":
		ksVal = config.get("lens_settings", "nfw_ks_val")
		intUpperLimit = config.get("lens_settings", "int_upper_limit") 
		prec = config.get("lens_settings", "prec")
		return([ksVal, intUpperLimit, prec])

def ampfachandler(config, injection_parameters, w_array_file, y_array_file, lens_model, additional_lens_parameters=[], mode="local"):
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
                amp_fac_real_file = data_subdir+"/fReal.dat"
                amp_fac_imag_file = data_subdir+"/fImag.dat" 

                if mode == "local":
                    print("Generating Lens Data") 
                    proc_to_run = [lens_model, w_array_file, y_array_file, amp_fac_real_file, amp_fac_imag_file] + additional_lens_parameters
                    subprocess.run(proc_to_run)
                    print("Lens Data Generated") 
                elif mode == "pipe":
                    print("Generating Lens Data Submit File")
                    generate_lens_subfile(config, amp_fac_real_file, amp_fac_imag_file, lens_model, w_array_file, y_array_file, additional_lens_parameters, outdir)
                    print("Lens Data Submit File Generated") 

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
        
        amp_fac_real = np.transpose(amp_fac_real)
        amp_fac_imag = np.transpose(amp_fac_imag)
        
        amp_interp_real = scint.RectBivariateSpline(w_array, y_array, amp_fac_real, kx=1, ky=1, s=0)
        amp_interp_imag = scint.RectBivariateSpline(w_array, y_array, amp_fac_imag, kx=1, ky=1, s=0)
        
        interpolator_func = lambda w, y: (amp_interp_real(w,y)+1j*amp_interp_imag(w,y)).flatten()
        
        print("Interpolator Generated")
        
        return(interpolator_func)

def generate_lens_subfile(config, amp_fac_real_file, amp_fac_imag_file, lens_model, w_array_file, y_array_file, additional_lens_parameters, outdir):
        submit_directory = outdir+"/submit" 
        if not os.path.isdir(submit_directory):
            os.mkdir(submit_directory) 
        subfile = submit_directory+"/generate_lens.sub"
        condor_settings_dictionary = config._sections["condor_settings"].copy()

        sub = open(subfile, "w")
        sub.write("universe = vanilla\n")
        sub.write("transfer_input_files = " + os.path.abspath(w_array_file) + "," + os.path.abspath(y_array_file)+"\n")
        sub.write("executable = "+config.get("lens_settings","executable_directory")+"/"+lens_model+"\n") 
        
        arguments = os.path.abspath(w_array_file) + " " + os.path.abspath(y_array_file) + " " + os.path.abspath(amp_fac_real_file) + " " + os.path.abspath(amp_fac_imag_file) 
        for i in additional_lens_parameters:
            arguments = arguments + " " + i

        sub.write("arguments = " + arguments+"\n")
        sub.write("log = " + submit_directory + "/lens_generation.log\n")
        sub.write("output = " + submit_directory + "/lens_generation.out\n")
        sub.write("error = " + submit_directory + "/lens_generation.err\n")

        sub.write("should_transfer_files = " + condor_settings_dictionary["transfer_files"] + "\n")
        sub.write("when_to_transfer_output = " + condor_settings_dictionary["when_to_transfer_output"] + "\n")
        
        checklist = ["request_cpus", "request_memory", "accounting_group", "accounting_group_user"]
        for item in checklist:
            if item in condor_settings_dictionary:
                sub.write(item+" = "+condor_settings_dictionary[item]+"\n")

        sub.write("queue 1\n") 

def gen_inject_file(config, outdir, injection_parameters):
        print("Generating Injection File")     

        inject_filename = outdir+"/data/injection.prior"
        inject_dat_filename = outdir+"/data/injection.dat" 
        inject_file = open(inject_filename, "w")

        for key, value in injection_parameters.items():
            inject_file.write(str(key) + "=" + str(value) + "\n")

        inject_file.close() 

        subprocess.run(["bilby_pipe_create_injection_file", inject_filename, "--n-injection", config.get("bilby_pipe_settings", "n-injections"), "-f", inject_dat_filename])

        print("Injection File Generated") 

        return(inject_dat_filename) 

def gen_bilby_pipe_ini(config, outdir, inject_file, waveform_arguments):
        bilby_pipe_ini_filename = config.get("bilby_setup", "label") + "_bilby_pipe.ini"
        bilby_pipe_ini = open(bilby_pipe_ini_filename, "w")

        bilby_pipe_config = dict() 

        bilby_pipe_config["label"] = config.get("bilby_setup","label")
        bilby_pipe_config["outdir"] = outdir
        
        detectors = config.get("bilby_setup","detectors")
        detector_list = list(detectors.split(","))
        bilby_pipe_config["detectors"] = detector_list 

        bilby_pipe_config["duration"] = config.get("bilby_setup", "duration")
        bilby_pipe_config["sampling_frequency"] = config.get("bilby_setup", "sampling_frequency")

        bilby_pipe_config["injection"] = config.get("bilby_pipe_settings", "injection")
        bilby_pipe_config["injection_file"] = inject_file

        bilby_pipe_config["gaussian-noise"] = config.get("bilby_pipe_settings", "gaussian-noise")
        bilby_pipe_config["waveform-generator"] = config.get("bilby_pipe_settings", "waveform-generator-class") 
        bilby_pipe_config["frequency-domain-source-model"] = config.get("bilby_pipe_settings", "frequency-domain-source-model")

        for key, value in config._sections["condor_settings"].items():
            if key != "when_to_transfer_output": 
                bilby_pipe_config[key] = value 

        bilby_pipe_config["sampler"] = config.get("sampler_settings", "sampler")

        sampler_kwargs_dict = dict() 
        for key, value in config._sections["sampler_settings"].items():
            if key != "sampler":
                sampler_kwargs_dict[key] = value
        bilby_pipe_config["sampler-kwargs"] = sampler_kwargs_dict 

        bilby_pipe_config["prior-file"] = config.get("prior_settings", "prior_file")

        bilby_pipe_config["waveform_arguments_dict"] = waveform_arguments

        bilby_pipe_configparser = ConfigParser()
        bilby_pipe_configparser["DEFAULT"] = bilby_pipe_config
        bilby_pipe_configparser.write(bilby_pipe_ini) 

        bilby_pipe_ini.close() 

        return()

def gen_overarch_dag(config, outdir):
        submit_directory = outdir+"/submit"
        lens_generation_subfile = submit_directory + "/generate_lens.sub"
        for filename in os.listdir(submit_directory):
            if filename.startswith("dag"):
                bilby_pipe_dag = submit_directory + "/" + filename

        overarch_filename = submit_directory + "/dag_" + config.get("bilby_setup", "label") + "_overarch.submit" 

        overarch_dag = open(overarch_filename, "w")

        overarch_dag.write("JOB lens_generation " + os.path.abspath(lens_generation_subfile) + "\n")
        overarch_dag.write("SUBDAG EXTERNAL bilby_pipe " + os.path.abspath(bilby_pipe_dag) + "\n")
        overarch_dag.write("PARENT lens_generation CHILD bilby_pipe")

        overarch_dag.close() 

        return() 
