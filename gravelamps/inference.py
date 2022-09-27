'''
Gravelamps Inference Run

Functions within control creating an inference run using bilby
The main function will get settings from the INI and will perform the run

Written by Mick Wright 2022
'''

from argparse import Namespace
import importlib
import os

from bilby_pipe.main import (MainInput,
                             perform_runtime_checks,
                             write_complete_config_file)
from bilby_pipe.job_creation import generate_dag
import htcondor

from gravelamps.core.gravelog import gravelogger, setup_logger
from gravelamps.core.graveparser import (create_graveparser,
                                         get_bilbypipe_args)
from gravelamps.core.file_handling import (create_bilby_pipe_config,
                                           create_final_dag,
                                           create_injection_file,
                                           get_config,
                                           get_output_directories)
from gravelamps.core.module_handling import get_lens_module
from gravelamps.generate_lens import main as generate_lens

def generate_files(config, lens_module, lens_package, lens_generation_args, injection):
    '''
    Input:
        lens_module - Name of module containing lensing functions
        lens_package - Module containing lensing functions
        config - INI configuration parser
        lens_generation_args - Commandline arguments passed to program
        injection - Boolean flag to indicate whether to use injection or analysis arguments

    Output:
        interpolator_files - Dictionary of files for interpolator modules, None otherwise

    Function generates files necessary for interpolator if the module requires it, returning None
    otherwise
    '''

    if lens_module == "gravelamps.lensing.interpolator":
        interpolator_files = generate_lens(_config=config,
                                           _args=lens_generation_args,
                                           _injection=injection)
    elif hasattr(lens_package, "generate_interpolator_files"):
        interpolator_files = lens_package.generate_interpolator_files(config,
                                                                      lens_generation_args,
                                                                      injection)
    else:
        interpolator_files = None

    return interpolator_files

def generate_waveform_arguments(config, args):
    '''
    Input:
        config - INI configuration parser
        args - Commandline arguments passed to the program

    Ouput:
        waveform_arguments - Dictionary of arguments for the waveform
    '''

    lens_module = get_lens_module(config, args)
    lens_package = importlib.import_module(lens_module)

    gravelogger.info("Lens Module: %s", lens_module)

    interpolator_files = generate_files(config, lens_module, lens_package,  args, args.injection)

    waveform_arguments = {}
    if interpolator_files is not None:
        for key, value in interpolator_files.items():
            if key == "source_position":
                waveform_arguments["source_position_file"] = os.path.abspath(value)
            else:
                waveform_arguments[key] = os.path.abspath(value)

    if args.injection:
        prepend = "injection_"
    else:
        prepend = ""

    parameter_list = ("waveform_approximant", "reference_frequency",
                      "minimum_frequency", "maximum_frequency")
    for parameter in parameter_list:
        value = config.get("inference_settings", f"{prepend}{parameter}", fallback=None)
        if value is None:
            value = config.get("inference_settings", f"{parameter}")
        waveform_arguments[parameter] = value

    if lens_module == "gravelamps.lensing.nfw":
        waveform_arguments["scaling_constant"] = config.get(f"{prepend}lens_generation_settings",
                                                            "nfw_scaling_constant",
                                                            fallback=None)
        if waveform_arguments["scaling_constant"] is None:
            waveform_arguments["scaling_constant"] = config.get("lens_generation_settings",
                                                                "nfw_scaling_constant")
    elif lens_module == "gravelamps.lensing.o3point":
        waveform_arguments["lookup_table_location"] =\
            config.get("run_settings", "lookup_table_location")

    elif lens_module == "gravelamps.lensing.millilensing":
        waveform_arguments["millilensing_kmax"] = config.get("run_settings", "millilensing_kmax")

    waveform_arguments["lens_module"] = lens_module

    return waveform_arguments

def generate_inference_args(args, injection=False):
    '''
    Input:
        args - Commandline arguments to the program
        injection - Boolean flag to set args.injection

    Output:
        inference_args - Modified version of submitted arguments

    Function creates a version of the arguments provided to inference to be used as part of the
    inference script
    '''

    inference_args = Namespace(**vars(args))

    inference_args.injection = injection
    inference_args.submit = False

    print(inference_args)

    return inference_args

def run_bilby_pipe_functions(bilby_pipe_config):
    '''
    Input:
        bilby_pipe_config - Dictionary containing options for the running of bilby_pipe functions

    Function runs bilby_pipe functions to generate data generation and analysis jobs.
    '''

    bilby_pipe_parser, bilby_pipe_args, bilby_pipe_unknown_args =\
        get_bilbypipe_args(bilby_pipe_config)
    bilby_pipe_input = MainInput(bilby_pipe_args, bilby_pipe_unknown_args)
    perform_runtime_checks(bilby_pipe_input, bilby_pipe_args)
    bilby_pipe_input.pretty_print_prior()
    write_complete_config_file(bilby_pipe_parser, bilby_pipe_args, bilby_pipe_input)
    generate_dag(bilby_pipe_input)

def main():
    '''
    Function will read in the commandline arguments and configuration INI and use this to create a
    lensed inference run using bilby.
    '''

    graveparser = create_graveparser()
    args = graveparser.parse_args()

    config = get_config(args)

    if config.getboolean("run_settings", "injection"):
        args.injection = True
    if config.getboolean("run_settings", "local"):
        args.local = True

    output_directories = get_output_directories(config)
    logging_level = config.get("output_settings", "logging_level", fallback="INFO")
    setup_logger(output_directories["outdir"], logging_level, args=args)
    gravelogger.info("Gravelogger setup complete, log will be output to %s/grave.log",
                     output_directories["outdir"])

    if args.injection or config.getboolean("run_settings", "injection"):
        gravelogger.info("Generating Injection Waveform Arguments")
        injection_args = generate_inference_args(args, injection=True)
        injection_waveform_arguments = generate_waveform_arguments(config, injection_args)
        gravelogger.info("Injection Waveform Arguments: %s", injection_waveform_arguments)
        gravelogger.info("Generating Injection")
        injection_file = create_injection_file(config)
        gravelogger.info("Injection File Generated")
    else:
        injection_file = None
        injection_waveform_arguments = None

    gravelogger.info("Generating Analysis Waveform Arguments")
    analysis_args = generate_inference_args(args, injection=False)
    analysis_waveform_arguments = generate_waveform_arguments(config, analysis_args)
    gravelogger.info("Analysis Waveform Argumnets: %s", analysis_waveform_arguments)

    gravelogger.info("Running bilby_pipe functions")
    bilby_pipe_config =\
        create_bilby_pipe_config(config, args, output_directories,
                                 injection_file=injection_file,
                                 injection_waveform_arguments=injection_waveform_arguments,
                                 analysis_waveform_arguments=analysis_waveform_arguments)
    run_bilby_pipe_functions(bilby_pipe_config)
    gravelogger.info("Bilby Pipe functions completed")

    if not args.local or config.getboolean("run_settings", "local", fallback=False):
        gravelogger.info("Constructing final DAG")
        dag_file = create_final_dag(config, output_directories)
        if args.submit or config.getboolean("run_settings", "submit", fallback=False):
            gravelogger.info("Submitting DAG to scheduler")
            htcondor.Submit.from_dag(dag_file)
        else:
            gravelogger.info("To submit, use the following command: \n\
                             \t$ condor_submit_dag %s", dag_file)

if __name__ == "__main__":
    main()
