'''
Gravelamps File Handling

Functions within handle input and output files from the Gravelamps programs.

Written by Mick Wright 2022
'''

from configparser import ConfigParser
import os
import shutil

import numpy as np

from bilby_pipe.create_injections import create_injection_file as bilby_injection_file
from bilby_pipe.input import Input

from gravelamps.core.gravelog import gravelogger

def create_bilby_pipe_config(config, args, output_directories, **kwargs):
    '''
    Input:
        config - INI configuration parser
        args - Commandline arguments given to program
        output_directories - Dictionary containing output directories for run

        **kwargs:
            injection_file - Location of file containing injection data
            analysis_waveform_arguments - Analysis waveform arguments
            injection_waveform_arguments - Injection waveform arguments

    Output:
        bilby_pipe_config - Dictionary containing settings for bilby_pipe

    Function generates a bilby_pipe configuration dictionary
    '''

    additional_information = {"injection_file": None,
                              "analysis_waveform_arguments": None,
                              "injection_waveform_arguments": None}

    additional_information.update(**kwargs)

    gravelogger.info("Generating bilby_pipe configuration")

    bilby_pipe_config = {}

    bilby_pipe_config["local"] = args.local or config.getboolean("run_settings",
                                                                 "local",
                                                                 fallback=False)
    bilby_pipe_config["submit"] = False

    bilby_pipe_config["accounting"] = config.get("condor_settings", "accounting_group", fallback=None)
    if bilby_pipe_config["accounting"] is None:
        gravelogger.warning("No accounting tag, jobs will not submit!")

    for key, value in config.items("condor_settings"):
        if key == "accounting":
            continue

        if "request_" in key:
            bilby_pipe_config[key] = value.replace("GB", "").replace("MB", "")
        else:
            bilby_pipe_config[key] = value

    bilby_pipe_config["label"] = config.get("output_settings", "label")
    bilby_pipe_config["outdir"] = output_directories["outdir"]

    for key, value in config.items("inference_settings"):
        if key == "prior":
            bilby_pipe_config["prior-file"] = value
        else:
            bilby_pipe_config[key.replace("_","-")] = value

    for key, value in config.items("bilby_pipe_additional_settings"):
        bilby_pipe_config[key.replace("_", "-")] = value

    if "overwrite-outdir" not in bilby_pipe_config:
        bilby_pipe_config["overwrite-outdir"] = True

    bilby_pipe_config["injection"] = args.injection or\
                                     config.getboolean("run_settings", "injection")

    if bilby_pipe_config["injection"]:
        bilby_pipe_config["gaussian-noise"] = True
        bilby_pipe_config["injection_file"] = additional_information["injection_file"]
        bilby_pipe_config["injection_waveform_arguments"] =\
            additional_information["injection_waveform_arguments"]

    bilby_pipe_config["waveform_arguments_dict"] =\
        additional_information["analysis_waveform_arguments"]

    gravelogger.info("bilby_pipe configuration generated")

    return bilby_pipe_config

def create_final_dag(config, output_directories):
    '''
    Input:
        config - INI configuration parser
        output_directories - Dictionary containing output directories

    Output:
        final_dag - Final DAG fila containing the linked jobs for inference runs

    Function takes the submit files created by Gravelamps as well as the DAG file made by
    bilby_pipe and links them together with the correct parent-child relations
    '''

    label = config.get("output_settings", "label")

    final_dag_file = f"{output_directories['submit']}/gravelamps_inference.dag"

    injection_generation_file = f"{output_directories['submit']}/generate_injection_interpolator_data.sub"
    analysis_generation_file = f"{output_directories['submit']}/generate_analysis_interpolator_data.sub"
    bilby_pipe_dag_file = f"{output_directories['submit']}/dag_{label}.submit"

    with open(final_dag_file, "w", encoding="utf-8") as dag:
        if os.path.isfile(injection_generation_file):
            gravelogger.info("Adding injection lens generation to final DAG")
            dag.write(f"JOB injection_lens_generation {injection_generation_file} \n")
        if os.path.isfile(analysis_generation_file):
            gravelogger.info("Adding analysis lens generation to final DAG")
            dag.write(f"JOB analysis_lens_generation {analysis_generation_file} \n")
        gravelogger.info("Adding bilby_pipe inference to final DAG")
        dag.write(f"SUBDAG EXTERNAL bilby_pipe_dag {bilby_pipe_dag_file} \n")

        if os.path.isfile(injection_generation_file):
            gravelogger.info("Parent-Child linking injection lens generation")
            dag.write("PARENT injection_lens_generation CHILD bilby_pipe_dag \n")
        if os.path.isfile(analysis_generation_file):
            gravelogger.info("Parent-Child linking analysis lens generation")
            dag.write("PARENT analysis_lens_generation CHILD bilby_pipe_dag \n")

    gravelogger.info("DAG file generated")

    return final_dag_file

def create_injection_file(config):
    '''
    Input:
        config - INI configuration parser
        injection_args - Arguments passed to program

    Output:
        injection_file - File containing injection data

    Function creates an injection file using bilby_pie_create_injection_file functions
    '''

    injection_file = f"{config.get('output_settings', 'outdir')}/data/injection.dat"

    prior_dict = config.items("injection_parameters")

    prior_file = injection_file.replace("injection", "prior")
    with open(prior_file, "w", encoding="utf-8") as prior:
        for key, value in prior_dict:
            prior.write(f"{key} = {value} \n")

    gpstuple = config.get("bilby_pipe_additional_settings", "gps_tuple", fallback=None)
    gpsfile = config.get("bilby_pipe_additional_settings", "gps_file", fallback=None)

    if gpstuple is not None:
        gpstimes = Input.parse_gps_tuple(gpstuple)
    elif gpsfile is not None:
        gpstimes = Input.read_gps_file(gpsfile)
    else:
        gpstimes=None

    n_injection = config.getint("bilby_pipe_additional_settings", "n-simulation")
    trigger_time = config.getfloat("inference_settings", "trigger_time", fallback=0)
    delta_t = config.getfloat("bilby_pipe_additional_settings", "deltaT", fallback=0.2)
    duration = config.getfloat("inference_settings", "duration", fallback=4)
    post_trigger_duration = config.getfloat("bilby_pipe_additional_settings",
                                            "post_trigger_duration",
                                            fallback=2)
    generation_seed = config.getfloat("bilby_pipe_additional_settings",
                                      "generation_seed",
                                      fallback=None)
    enforce_signal_duration = config.getboolean("bilby_pipe_additional_settings",
                                                "enforce_signal_duration",
                                                fallback=False)

    bilby_injection_file(injection_file,
                         prior_file=prior_file,
                         prior_dict=None,
                         n_injection=n_injection,
                         trigger_time=trigger_time,
                         deltaT=delta_t,
                         gpstimes=gpstimes,
                         duration=duration,
                         post_trigger_duration=post_trigger_duration,
                         generation_seed=generation_seed,
                         enforce_signal_duration=enforce_signal_duration)

    return injection_file

def get_config(args):
    '''
    Input:
        args - Commandline arguments given to program

    Output:
        config - INI configuration parser

    Function retrieves the user defined INI file and loads it into a ConfigParser object.
    Will raise errors if the user has not specified the file or the file cannot be read.
    '''

    config = ConfigParser()

    if not os.path.isfile(args.ini):
        raise IOError(f"{args.ini} is not a valid file")

    try:
        config.read(args.ini)
    except IOError:
        print(f"{args.ini} cannot be read!")

    return config

def get_output_directories(config, from_config=True):
    '''
    Input:
        config - INI configuration parser
        from_config - Boolean value to indicate whether to use config value or assume being run
                      the data subdirectory

    Output:
        output_dir_dict - Dictionary containing output directories for run

    Function gets the output directory from the INI and the relevant subdirectories for the run
    type i.e. local or not. Function will create directory structure if it does not already exist.
    '''

    if from_config:
        outdir = config.get("output_settings", "outdir", fallback=".")
    else:
        outdir = os.path.dirname(os.getcwd())

    data_subdirectory = f"{outdir}/data"
    submit_subdirectory = f"{outdir}/submit"

    output_dir_dict = {"outdir": outdir, "data": data_subdirectory, "submit": submit_subdirectory}

    for _, value in output_dir_dict.items():
        if not os.path.isdir(value):
            os.mkdir(value)

    return output_dir_dict

def retrieve_interpolator_files(config, args):
    '''
    Input:
        config - INI configuration parser
        args - Commandline arguments given to program

    Output:
        file_dict - dictionary of files containing dimesnionlesss frequency and source position
                    values over which to generate the interpolator, followed by the corresponding
                    files containin the real and imaginary amplification factor values to use as
                    the interpolating data

    Function retrieves from the INI the locations of the files containing the data necessary to
    generate the interpolator, checking that these files exist throwing errors if not. Will proceed
    normally if files are missing, leaving them omitted from the resulting dictionary.
    '''

    if args.injection:
        lens_type = "injection"
    else:
        lens_type = "analysis"

    value_types = ("dimensionless_frequency", "source_position",
                   "amplification_factor_real", "amplification_factor_imag")

    file_dict = {}
    for value in value_types:
        file_location = config.get(f"{lens_type}_lens_generation_settings",
                                   f"{value}_file", fallback=None)
        if file_location is None or file_location == "None":
            gravelogger.info("%s file not given", value.replace("_", " "))
            file_location = None
        else:
            if not os.path.isfile(file_location):
                raise IOError(f"{file_location} not found!")
            gravelogger.info("%s file found at %s", value.replace("_", " ").title(), file_location)
        file_dict[value] = file_location

    gravelogger.info("Files found before handling: %s", file_dict)
    return file_dict

def grid_file_handler(config, args, data_subdirectory, file_dict):
    '''
    Input:
        config - INI configuration parser
        args - Commandline arguments given to program
        data_subdirectory - Subdirectory of run output directory containing data files
        file_dict - Dictionary of files containing either location of files containing the
                    dimensionless frequency and source position or Nones to indicate these
                    files require generation

    Output:
        lens_file_dict - Dictionary containing both grid files needed for the interpolator
                         may be equivalent to file_dict if files exist already and the don't
                         make copy option is selected in config

    Function handles the dimensionless frequency and source position value files - i.e. the grid
    axis files for the interpolators. If the user has specified them and to copy, it will copy these
    files to the default locations. It will then return the final location of the files in a
    dictionary.
    '''

    if args.injection:
        lens_type = "injection"
    else:
        lens_type = "analysis"

    amplification_defined = bool(file_dict["amplification_factor_real"])\
                            or bool(file_dict["amplification_factor_imag"])

    temp_dict = {}
    for file_type in ("dimensionless_frequency", "source_position"):
        default_outfile = f"{data_subdirectory}/{lens_type}_{file_type}.dat"

        if file_dict[file_type] is None:
            if amplification_defined:
                gravelogger.warning(("Amplification factor files defined without corresponding "\
                                     "%s file. Interpolator grid may not be accurate!"),
                                     file_type.replace("_", " "))

            gravelogger.info("Generating %s file", file_type.replace("_", " "))

            min_value = config.getfloat(f"{lens_type}_lens_generation_settings",
                                        f"minimum_{file_type}")
            max_value = config.getfloat(f"{lens_type}_lens_generation_settings",
                                        f"maximum_{file_type}")
            num_values = config.getint(f"{lens_type}_lens_generation_settings",
                                       f"length_{file_type}")

            value_array = np.linspace(min_value, max_value, num_values)
            np.savetxt(default_outfile, value_array)
            temp_dict[file_type] = default_outfile

        else:
            shutil.copyfile(file_dict[file_type], default_outfile)
            temp_dict[file_type] = default_outfile

    lens_file_dict = dict(file_dict, **temp_dict)
    gravelogger.info("Files after grid handling: %s", lens_file_dict)

    return lens_file_dict

def data_file_handler(args,
                      data_subdirectory,
                      file_dict):
    '''
    Input:
        args - Commandline arguments given to program
        data_subdirectory - Subdirectory of run output directory containing data files
        file_dict - Dictionary containing the list of file locations for constructing
                    the interpolator before data file handling

    Output:
        lens_file_dict - Dictionary containing the list of file locations for constructing
                         the interpolator after data file handling
        complete_files - Number of the data files that are complete
    '''

    if args.injection:
        lens_type = "injection"
    else:
        lens_type = "analysis"

    dimensionless_frequency_array = np.loadtxt(file_dict["dimensionless_frequency"])
    source_position_array = np.loadtxt(file_dict["source_position"])

    grid_shape = (len(dimensionless_frequency_array), len(source_position_array))

    temp_dict = {}
    complete_files = 0
    for file_type in "amplification_factor_real", "amplification_factor_imag":
        default_outfile = f"{data_subdirectory}/{lens_type}_{file_type}.dat"

        if bool(file_dict[file_type]):
            data_grid = np.loadtxt(file_dict[file_type])

            if grid_shape in (data_grid.shape, data_grid.transpose.shape):
                complete_files += 1

            shutil.copyfile(file_dict[file_type], default_outfile)
            temp_dict[file_type] = default_outfile
        else:
            temp_dict[file_type] = default_outfile

    lens_file_dict = dict(file_dict, **temp_dict)

    gravelogger.info("Files after data file handling: %s", lens_file_dict)
    gravelogger.info("%s data file(s) complete", complete_files)

    return lens_file_dict, complete_files
