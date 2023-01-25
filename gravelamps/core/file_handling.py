"""Gravelamps File Handling

Following are functions handling input and output files from the Gravelamps programs.

Written by Mick Wright 2022

Routines
--------
create_bilby_pipe_config
    Generate a bilby_pipe configuration dictionary
create_final_dag
    Generate the overall gravelamps DAG
create_injection_file
    Generate a bilby_pipe injection file
get_config
    Retrieves user INI configuration from arguments
get_output_directories
    Retrieves the output directories
retrieve_interpolator_files
    Retrieves the files necessary for lens interpolator generation
grid_file_handler
    Handles the interpolator grid file generation and locations
data_file_handler
    Handles the interpolator data file generation and locations

"""

from configparser import ConfigParser
import os
import shutil

import numpy as np

from bilby_pipe.create_injections import create_injection_file as bilby_injection_file
from bilby_pipe.input import Input

from gravelamps.core.gravelog import gravelogger

def create_bilby_pipe_config(config, args, output_directories, **kwargs):
    """Generates a bilby_pipe configuration dictionary

    Parameters
    ----------
    config : configparser.ConfigParser
        Object containing settings from INI file
    args : argparse.Namespace
        Object containing commandline arguments to the program
    output_directories : dict
        Contains the output directories for the run
    injection_file : str, optional
        Path of file containing injection data
    analysis_waveform_arguments : dict, optional
        Arguments dictionary to the analysis waveform generator
    injection_waveform_arguments : dict, optional
        Arguments dictionary for the injection waveform generator

    Returns
    -------
    bilby_pipe_config : dict
        Contains the configuration settings for a bilby_pipe run
    """

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

    bilby_pipe_config["accounting"] =\
        config.get("condor_settings", "accounting_group", fallback=None)
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
    """
    Generate overall Gravelamps DAG.

    This DAG will contain the jobs to be submitted to the HTCondor scheduler with correct parent
    child linking. Lens generation jobs will be run first and may run with no linking to each
    other. These jobs form the parents of the bilby_pipe inference runs using the lensed waveforms.

    Parameters
    ----------
    config : configparser.ConfigParser
        Object containing settings from INI file
    output_directories : dict
        Contains the output directories for the run

    Returns
    -------
    final_dag : str
        Path to the gravelamps DAG file
    """

    label = config.get("output_settings", "label")

    final_dag_file = f"{output_directories['submit']}/gravelamps_inference.dag"

    injection_generation_file =\
        f"{output_directories['submit']}/generate_injection_interpolator_data.sub"
    analysis_generation_file =\
        f"{output_directories['submit']}/generate_analysis_interpolator_data.sub"
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
    """
    Generate bilby_pipe injection file.

    Parameters
    ----------
    config : configparser.ConfigParser
        Object containing settings from INI file

    Returns
    -------
    injection_file : str
        Path to the created bilby_pipe injection file
    """

    injection_file = f"{config.get('output_settings', 'outdir')}/data/injection.dat"
    prior_dict = config.items("injection_parameters")
    prior_file = injection_file.replace("injection.dat", "prior.dat")

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
    """
    Retrieves user INI configuration from arguments

    Parameters
    ----------
    args : argparse.Namespace
        Object containing commandline arguments to program

    Returns
    -------
    config : configparser.ConfigParser
        Object containing settings from INI file

    Raises
    ------
    IOError
        Where the INI file is not specified within the arguments or cannot be read
    """

    config = ConfigParser()

    if not os.path.isfile(args.ini):
        raise IOError(f"{args.ini} is not a valid file")

    try:
        config.read(args.ini)
    except IOError:
        print(f"{args.ini} cannot be read!")

    return config

def get_output_directories(config, from_config=True):
    """
    Retrieves the output directories.

    The output directories specified are the top level output directory, followed by data and submit
    subdirectories with the specified names, 'data', and 'submit'.

    The top level directory is typically specified within the user specified INI file. This will
    fallback to the current working directory, or can be specified to directly run presuming such.
    These folders will be created if they are not already extant.

    Parameters
    ----------
    config : configparser.ConfigParser
        Object containing settings from INI file
    from_config : bool, optional
        Flag to ignore the INI and set the top level directory to the current directory

    Returns
    -------
    output_dir_dict : dict
        Contains the `output` top level directory and `submit` and `data` subdirectories in
        the specified keys.
    """

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
    """
    Retrieves files necessary for generation of the lens interpolator

    Will proceed if the file does not exist---specifying that it needs to be created. Will throw
    exception if a file that does not exist is specified to exist.

    Parameters
    ----------
    config : configparser.ConfigParser
        Object containing settings from INI file
    args : argparse.Namespace
        Object containing commandline arguments to program

    Returns
    -------
    file_dict : dict
        Contains `dimnesionless_frequency`, `source_position`, `amplification_factor_real`,
        `amplification_factor_imag` keys. Each of these is a string path to file containing the
        specified data for the lens interpolator generation.

    Raises
    ------
    IOError
        In case where a file is specified in the INI that does not exist
    """

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
    """
    Handles the interpolator grid file generation and locations

    These files specify the dimensionless frequency and source position grid structure that is
    interpolated over for the amplification factor data. These files may be directly specified
    in the INI, or may be constructed if not. User will be warned if the amplification factor
    files are defined without also defining the grid files, since the grid may not be accurate
    if generated for pre-existing data.

    Parameters
    ----------
    config : configparser.ConfigParser
        Object containing settings from INI file
    args : argparse.Namesapce
        Òbject containing commandline arguments to program
    data_subdirectory : str
        Path to the subdirectory containing data files
    file_dict : dict
        Contains either location of grid files, or None indicating these files require generation

    Returns
    -------
    lens_file_dict : dict
        Contains locations of completed grid files. Equivalent to file_dict if these are specified
        in the data_subdirectory.
    """

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
    """
    Handles the interpolator data file generation and locations.

    These files specify the real and imaginary components of the amplification factor data that
    forms the base of the interpolator objects. These files may be directly specified in the INI
    or may not be, specifying that they need generation. This handler does not run the generation
    itself due to the computational complexity, instead it specifies the number of files that are
    complete.

    Parameters
    ----------
    args : argparse.Namespace
        Object containing commandline arguments to program
    data_subdirectory : str
        Path to the subdirectory containing the data files
    file_dict : dict
        Contains either the path to the files, or None to indicate they require generation

    Returns
    -------
    lens_file_dict : dict
        Contains the path to the files for construction of lens interpolator. Equivalent to
        `file_dict` if these files are specified in `data_subdirectory`
    complete_files : int
        Number of data files that are complete
    """

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

            if grid_shape in (data_grid.shape, data_grid.transpose().shape):
                complete_files += 1

            shutil.copyfile(file_dict[file_type], default_outfile)
            temp_dict[file_type] = default_outfile
        else:
            temp_dict[file_type] = default_outfile

    lens_file_dict = dict(file_dict, **temp_dict)

    gravelogger.info("Files after data file handling: %s", lens_file_dict)
    gravelogger.info("%s data file(s) complete", complete_files)

    return lens_file_dict, complete_files
