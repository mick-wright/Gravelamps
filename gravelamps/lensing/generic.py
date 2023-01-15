'''
Generic Lensing Functions

These functions are generic non-model-specific functions that assist in the running of the model
specific code.

Written by Mick Wright 2022
'''

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
    '''
    Input:
        config - INI configuration parser
        args - Commandline arguments passed to the program
        argument_list - list of the model specific arguments to retrieve from  the config file
        type_list - list of the types of the arguments above

    Output:
        value_list - list of the argument values as the types specified
    '''

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
    '''
    Input:
        config - INI configuration parser
        args - Commandline arguments passed to the program
        output_directories - dictionary containing paths to the folders for output
        model - String containing the python path to the model
        file_dict - dictionary of files containing the paths for the interpolator grid and data

    Output:
        condor_settings - Dictionary containing condor settings to use

    Function generates the settings for the HTCondor Job based upon a set of defaults replaced by
    any user specified settings.
    '''

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
    '''
    Input:
        args - Commandline arguments passed to the program
        output_directories - dictionary containing paths to the folders for output
        file_dict - Files containing interpolator data grid and value files
        condor_settings - dictionary containing settings for the HTCondor job
    '''

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
    '''
    Input:
        config - INI configuration parser
        args - Commandline arguments passed to the program
        model - String containing the python path to the lensing module
        file_dict - dictionary of files containing the paths for the interpolator grid and data

    Function generates interpolator data using the specified model and saving the data to the given
    files.
    '''
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
    '''
    Runs generate_interolator_data with the specified model and file dictionary
    '''

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
