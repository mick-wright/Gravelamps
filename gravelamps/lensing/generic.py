"""Generic Lensing Functions

Following are functions that are model agnostic and used to interface the specific model forms with
the rest of the architecture. Also included is a main running script which forms the minimal program
`gravelamps_generate_interpolator_data` which is included for HTCondor scheduling of lens data
generation. This is not generally intended to be run by the user.

Written by Mick Wright 2022

Routines
--------
get_additional_arguments
    Retrieves required additional arguments for model specific lens generation
get_condor_config
    Generates condor settings from the user INI
generate_interpolator_condor
    Builds and, if specified, submits a lens generation run for the HTCondor scheduler
generate_inerpolator_data
    Performs lens generation for the specified model in the format needed
"""

import ast
import importlib
import os
import sys

import htcondor
import numpy as np

from gravelamps.core.file_handling import get_config, get_output_directories
from gravelamps.core.gravelog import gravelogger, setup_logger
from gravelamps.core.graveparser import create_graveparser

def get_additional_arguments(config, args, argument_list, type_list):
    """
    Retrieves required additional arguments for model specific lens generation

    Parameters
    ----------
    config : configparser.ConfigParser
        Object containing settings for user INI file
    args : argparse.Namespace
        Object containing commmandline arguments to program
    argument_list : List of str
        Model specific arguments to retrieve from config file
    type_list : List of type funcs
        Types of the arguments given in `argument_list`

    Returns
    -------
    value_list : List of objects
        Values of the arguments specified in `argument_list`
    """

    gravelogger.info("Retrieving additional arguments")

    if args.injection:
        lens_type = "injection"
    else:
        lens_type = "analysis"

    value_list = []
    for argument, typing in zip(argument_list, type_list):
        try:
            value = typing(config.get(f"{lens_type}_lens_generation_settings", argument))
        except ValueError:
            print(f"{argument} inexpressible as {typing}")
        value_list.append(value)

        gravelogger.info("Setting %s = %s", argument, value)

    return value_list

def get_condor_config(config, args, output_directories, model, file_dict):
    """
    Generates HTCondor configuration from user settings in INI

    Parameters
    ----------
    config : configparser.ConfigParser
        Object containing settings from user INI file
    args : argparse.Namespace
        Object containing commandline arguments to program
    output_directories : dict
        Contains paths to the folders for output
    model : str
        Full python path to the lens model
    file_dict : dict
        Contains the paths for the interpolator grid and data files

    Returns
    -------
    condor_settings : dict
        Contains settings used by HTCondor
    """

    if args.injection:
        lens_type = "injection"
    else:
        lens_type = "analysis"

    default_condor_settings = {
        "executable": f"{os.path.dirname(sys.executable)}/gravelamps_generate_interpolator_data",
        "initialdir": os.path.abspath(output_directories["data"]),
        "should_transfer_files": "YES",
        "when_to_transfer_output": "ON_EXIT_OR_EVICT",
        "output": f"{lens_type}_data_generation.out",
        "error": f"{lens_type}_data_generation.err",
        "log": f"{lens_type}_data_generation.log",
        "request_cpus": "16",
        "request_memory": "8 GB",
        "request_disk": "2 GB",
        }

    base_dict = {}

    for key, value in file_dict.items():
        if os.path.dirname(value) == output_directories["data"]:
            base_dict[key] = os.path.basename(value)
        else:
            base_dict[key] = os.path.abspath(value)

    dict_string = str(base_dict).replace("'", "####").replace(" ", "")
    default_condor_settings["arguments"] = f"{model} {dict_string} {os.path.abspath(args.ini)}"

    default_condor_settings["transfer_input_files"] =\
        f"{base_dict['dimensionless_frequency']}, {base_dict['source_position']},"\
	f"{base_dict['amplification_factor_real']}, {base_dict['amplification_factor_imag']}"

    user_condor_settings = dict(config.items("condor_settings"))
    if "accounting_group" not in user_condor_settings:
        gravelogger.warning("No accounting tag provided, condor submission may fail!")

    condor_settings = dict(default_condor_settings, **user_condor_settings)

    return condor_settings

def generate_interpolator_condor(args, output_directories, condor_settings):
    """
    Builds and, if specified, submits a lens generation run for the HTCondor Scheduler.

    Parameters
    ----------
    args : argparse.Namespace
        Object containing commandline arguments to program
    output_directories : dict
        Contains paths to the folders for output
    file_dict : dict
        Contains paths for interpolator grid and data files
    condor_settings : dict
        Contains the settings for HTCondor to use
    """

    if args.injection:
        lens_type = "injection"
    else:
        lens_type = "analysis"

    gravelogger.info("Generating %s interpolator data generation submission file", lens_type)
    submit_filename = f"{output_directories['submit']}/generate_{lens_type}_interpolator_data.sub"
    gravelogger.info("Submit file saved to %s", submit_filename)

    with open(submit_filename, "w", encoding="utf-8") as file:
        for key, value in condor_settings.items():
            file.write(f"{key} = {value} \n")
        file.write("queue 1")

    if args.submit:
        gravelogger.info("Submitting Job to Scheduler")

        condor_submit_object = htcondor.Submit(condor_settings)
        schedd = htcondor.Schedd()
        schedd.submit(condor_submit_object)

def generate_interpolator_data(config, args, model, file_dict):
    """
    Performs lens generation for the specified model in the format needed.

    Lens generation can be done in one of three ways --- it can be done locally from the
    `generate_interpolator_data` function for the model if it has it, it can be done locally
    from the `amplification_factor` function directly if the model is sufficiently speedy to not
    need a specific interpolator generation function, or it can generate an HTCondor DAG file to
    run one of these should scheduling be required.

    Parameters
    ----------
    config : configparser.ConfigParser
        Object containing settings from user INI file
    args : argparse.Namespace
        Object containing commandline arguments to program
    model : str
        Full python path to the lensing model
    file_dict : dict
        Contains paths for interpolator grid and data files

    Raises
    ------
    AttributeError
        Occurs when no data generation functions exist within model module
    """
    interpolator_module = importlib.import_module(model)

    if hasattr(interpolator_module, "generate_interpolator_data"):
        interpolator_module.generate_interpolator_data(config, args, file_dict)
    elif hasattr(interpolator_module, "amplification_factor"):
        gravelogger.info(("No specific interpolator generator function, using amplification_factor"\
                          " to generate data"))

        dimensionless_frequency_array = np.loadtxt(file_dict["dimensionless_frequency"])
        source_position_array = np.loadtxt(file_dict["source_position"])

        amplification_array =\
            interpolator_module.amplification_factor(dimensionless_frequency_array,
                                                     source_position_array)

        real_array = np.real(amplification_array)
        imag_array = np.imag(amplification_array)

        np.savetxt(file_dict["amplification_factor_real"], real_array)
        np.savetxt(file_dict["amplification_factor_imag"], imag_array)
    else:
        raise AttributeError(f"No data generating functions found in {model}")

def main():
    """
    Forms the program `gravelamps_generate_interpolator_data`.

    This program is not intended to be run directly by the user, instead it is a minimal function
    that will be called by the HTCondor scheduler when it is being used. Users should use the more
    full features `gravelamps_generate_lens` program. If run, however, it will use the specified
    file dictionary and output to generate lens data directly.
    """

    graveparser = create_graveparser()
    args = graveparser.parse_args()

    config = get_config(args)

    logging_level = config.get("output_settings", "logging_level", fallback="INFO")
    setup_logger(".", logging_level, args=args)
    gravelogger.info("Gravelogger setup complete, log will be output to grave.log")

    dictstring = args.filedict.replace("'","")
    dictstring = args.filedict.replace("####","'")
    dictstring = str(dictstring)
    file_dict = ast.literal_eval(dictstring)

    print(file_dict)

    generate_interpolator_data(config, args, args.model, file_dict)

if __name__ == "__main__":
    main()
