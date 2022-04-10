'''
Pipe file generators for gravelamps.inference

Functions generate the files for the pipe based inference.

Written by Mick Wright 2021
'''

import os
import subprocess

import bilby

def lens_subfile(config, dim_freq_file, sour_pos_file, amp_fac_real_file,
                 amp_fac_imag_file, additional_lens_parameters):
    '''
    Input:
        config - INI configuration parser
        dim_freq_file - location of file containing dimensionless frequency values
        sour_pos_file - location of file containing source position values
        amp_fac_real_file - location of file containing the real part of the amplification
                            factor values
        amp_fac_imag_file - location of file containing the imaginary part of the amplification
                            factor values
        additional_lens_parameters - list containing the additional parameters for the lens model

    Function generates the condor submission file in order to generate the amplification factor
    data necessary for the analysis run
    '''

    #Logging
    bilby.core.utils.logger.info("Generating Lens Generation Submission File") 

    #Get the submission and data subdirectories
    outdir = os.path.abspath(config.get("output_settings", "outdir"))
    submit_directory = f"{outdir}/submit"
    data_subdirectory = f"{outdir}/data"

    #Get the lens model
    lens_model = config.get("lens_generation_settings", "lens_model")

    #Get the executable directory and construct the full path of the lens model function
    executable_directory = config.get(
        "lens_generation_settings", "executable_directory", fallback="default")
    if executable_directory == "default":
        executable_directory = f"{os.path.expanduser('~')}/bin"
    executable = f"{executable_directory}/{lens_model}"

    #Create the submission subdirectory if it doesn't exist
    if not os.path.isdir(submit_directory):
        os.mkdir(submit_directory)

    #Name of submit file
    submit_file = f"{submit_directory}/lens_generation.sub"

    #Open the submit file and write the necessary parts
    with open(submit_file, "w", encoding="utf-8") as sub:
        sub.write("universe = vanilla\n")
        sub.write(f"initialdir = {data_subdirectory}\n")
        sub.write(f"executable = {executable}\n")

        #Construct the arguments necessary for the function call
        argument_list = [dim_freq_file, sour_pos_file, amp_fac_real_file, amp_fac_imag_file]\
                        + additional_lens_parameters

        sub.write(f"arguments = {' '.join(argument_list)}\n")
        sub.write(f"log = {data_subdirectory}/lens_generation.log\n")
        sub.write(f"output = {data_subdirectory}/lens_generation.out\n")
        sub.write(f"error = {data_subdirectory}/lens_generation.err\n")

        for key, value in config.items("condor_settings"):
            sub.write(f"{key} = {value}\n")

        sub.write("queue 1\n")

    #Logging End
    bilby.core.utils.logger.info("Lens Generation Submission File Generated") 

def injection_file(config, injection_parameters):
    '''
    Input:
        config - INI configuration parser
        injection_parameters = dictionary containing the injection parameters

    Output:
        inject_file - location of file containing the injection

    Function generates the necessary injection files for perform bilby_pipe analysis on
    injected simulated signals
    '''

    #Logging
    bilby.core.utils.logger.info("Generating Injection File")

    #Get the proper location for the injection file
    outdir = config.get("output_settings", "outdir")
    inject_prior = f"{outdir}/data/injection.prior"
    inject_file = f"{outdir}/data/injection.dat"

    #Write the prior file from the injection parameters
    with open(inject_prior, "w", encoding="utf-8") as prior:
        for key, value in injection_parameters.items():
            prior.write(f"{key} = {value}\n")

    #Run the bilby_pipe_create_injection_file function
    n_injections = config.get("injection_settings", "n-injection")

    subprocess.run(["bilby_pipe_create_injection_file", inject_prior,
                    "--n-injection", n_injections,
                    "-f", inject_file], check=True)

    #Logging End
    bilby.core.utils.logger.info("Injection File Generated")

    return inject_file

def bilby_pipe_ini(config, inject_file, injection_waveform_arguments, waveform_arguments, mode):
    '''
    Input:
        config - INI configuration parser
        injection_file - location of file containing injection data
        injection_waveform_arguments - dictionary containing arguments to be passed to waveform
                                       during data generatioon only
        waveform_arguments - dictionary containing arguments to be passed to waveform
        mode - string determinining whether to generate an unlensed or a lensed analysis run

    Output:
        ini_file - location of the generated INI file

    Function creates a bilby_pipe compatible INI from the options given by the user in the
    gravelamps INI
    '''

    #Get the bilby_pipe INI filename
    ini_file = config.get("output_settings", "label") + "_" + mode + "_bilby_pipe.ini"

    #Create the empty configuration dictionary
    bilby_pipe_config = dict()

    #Insert the accounting tag
    bilby_pipe_config["accounting"] = config.get("condor_settings", "accounting_group")

    #Read in the condor settings
    for key, value in config._sections["condor_settings"].items():
        bilby_pipe_config[key] = value

    #Insert the label and output directory
    bilby_pipe_config["label"] = config.get("output_settings", "label") + "_" + mode
    bilby_pipe_config["outdir"] = config.get("output_settings", "outdir")

    #Create the detector list
    detector_string = config.get("analysis_settings", "interferometers")
    detector_list = list(detector_string.split(","))

    bilby_pipe_config["detectors"] = detector_list

    #Insert the Waveform Generator Class and Frequency Domain Source Model depending on the mode
    if mode == "lensed":
        bilby_pipe_config["waveform-generator"] = config.get(
            "analysis_settings", "lensed_waveform_generator_class")
        bilby_pipe_config["frequency-domain-source-model"] = config.get(
            "analysis_settings", "lensed_frequency_domain_source_model")
    elif mode == "unlensed":
        bilby_pipe_config["waveform-generator"] = config.get(
            "unlensed_analysis_settings", "unlensed_waveform_generator_class")
        bilby_pipe_config["frequency-domain-source-model"] = config.get(
            "unlensed_analysis_settings", "unlensed_frequency_domain_source_model")

    #Include the duration and the sampling frequency
    bilby_pipe_config["duration"] = config.get("analysis_settings", "duration")
    bilby_pipe_config["sampling_frequency"] = config.get("analysis_settings", "sampling_frequency")

    #Include settings from bilby_pipe_settings
    for key, value in config._sections["bilby_pipe_settings"].items():
        bilby_pipe_config[key] = value

    #If injecting include all the necessary injection settings, otherwise include event settings
    if inject_file is not None:
        for key, value in config._sections["injection_settings"].items():
            bilby_pipe_config[key] = value
        bilby_pipe_config["injection_file"] = inject_file
        bilby_pipe_config["injection_waveform_arguments"] = injection_waveform_arguments 
    else:
        for key, value in config._sections["event_settings"].items():
            bilby_pipe_config[key] = value

    #Include the sampler settings
    bilby_pipe_config["sampler"] = config.get("analysis_settings", "sampler")
    bilby_pipe_config["sampler-kwargs"] = config._sections["sampler_kwargs"].copy()

    #Include the prior file
    if mode == "lensed":
        bilby_pipe_config["prior-file"] = config.get("analysis_settings", "prior_file")
    elif mode == "unlensed":
        bilby_pipe_config["prior-file"] = config.get("unlensed_analysis_settings", "prior_file")

    #Include the waveform arguments dictionary
    bilby_pipe_config["waveform_arguments_dict"] = waveform_arguments

    #Include the approximant specifically
    bilby_pipe_config["waveform_approximant"] = waveform_arguments["waveform_approximant"]

    #Include the trigger time
    bilby_pipe_config["trigger-time"] = config.get("analysis_settings", "trigger_time")

    #Include whether or not to plot the output as a corner plot
    bilby_pipe_config["plot-corner"] = config.get("analysis_settings", "plot_corner")

    #Write the dictionary to the file
    with open(ini_file, "w") as ini:
        for key, value in bilby_pipe_config.items():
            ini.write(key + " = " + str(value) + "\n")

    return ini_file

def overarching_dag(config):
    '''
    Input:
        config - INI configuration parser

    Output:
        dag_file - location of the file containing the overarching DAG that will submit the
                   lens generation, and bilby_pipe analysis runs

    Function generates the overarching DAG file to submit all of the condor jobs - the lensing
    gneration and the bilby_pipe analysis runs
    '''

    #Get the submission and data subdirectories
    outdir = config.get("output_settings", "outdir")
    label = config.get("output_settings", "label")
    submit_subdirectory = outdir + "/submit"

    #Get the Overarching DAG file name
    dag_file = submit_subdirectory + "/dag_" + label + "_overarch.submit"

    #Get the lens generation submit file
    lens_generation_subfile = submit_subdirectory + "/generate_lens.sub"
    lens_generation_subfile = os.path.abspath(lens_generation_subfile)

    #Open the DAG file for writing
    with open(dag_file, "w") as dag:
        #Determine methodology
        methodology = config.get("analysis_settings", "methodology")

        #If methodology is interpolate, and the lens files do not exist add the lens generation job
        if methodology == "interpolate":
            amp_fac_real_file = config.get(
                "lens_generation_settings", "amplification_factor_real_file")
            amp_fac_imag_file = config.get(
                "lens_generation_settings", "amplification_factor_imag_file")

            if not os.path.isfile(amp_fac_real_file) and not os.path.isfile(amp_fac_imag_file):
                dag.write("JOB lens_generation " + lens_generation_subfile + "\n")
                lens_generation = True
            else:
                lens_generation = False
        else:
            lens_generation = False

        #Add unlensed analysis run, if doing so
        if config.getboolean("unlensed_analysis_settings", "unlensed_analysis_run"):
            unlensed_sub_file = submit_subdirectory + "/dag_" + label + "_unlensed.submit"
            unlensed_sub_file = os.path.abspath(unlensed_sub_file)
            dag.write("SUBDAG EXTERNAL bilby_pipe_unlensed " + unlensed_sub_file + "\n")

        #Add lensed analysis run
        lensed_sub_file = submit_subdirectory + "/dag_" + label + "_lensed.submit"
        lensed_sub_file = os.path.abspath(lensed_sub_file)
        dag.write("SUBDAG EXTERNAL bilby_pipe_lensed " + lensed_sub_file + "\n")

        #Parent-Child Link the jobs so that the lens generation goes first, then the unlensed,
        #then the lensed analysis runs
        if config.getboolean("unlensed_analysis_settings", "unlensed_analysis_run"):
            if lens_generation is True:
                dag.write("PARENT lens_generation CHILD bilby_pipe_unlensed \n")
            dag.write("PARENT bilby_pipe_unlensed CHILD bilby_pipe_lensed \n")
        else:
            if lens_generation is True:
                dag.write("PARENT lens_generation CHILD bilby_pipe_lensed \n")

    return dag_file
