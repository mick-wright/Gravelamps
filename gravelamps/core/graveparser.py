'''
Gravelamps Argument Parser

Functions within control the use of the parser used within Gravelamps utilising the python bulitin
argparse module

Written by Mick Wright 2022
'''

import argparse
import json
import os
import sys

from bilby_pipe.parser import create_parser

def create_bilbypipe_parser():
    '''
    Output:
        bilbypipe_parser - Parser to be passed to bilby_pipe functions

    Function creates a parser for bilby_pipe's functions
    '''

    bilbypipe_parser = create_parser()

    return bilbypipe_parser

def create_graveparser():
    '''
    Output:
        graveparser - Gravelamps argument parser

    Creates argument parser based on the program being run
    '''

    if "generate_lens" in sys.argv[0]:
        graveparser = graveparser_lens_generation_parser()

    elif any(check in sys.argv[0]\
             for check in ("gravelamps_generate_interpolator_data", "generic.py", "condor_exec")):
        graveparser = graveparser_interpolator_parser()

    elif "inference" in sys.argv[0]:
        graveparser = graveparser_inference_parser()

    else:
        raise ValueError(f"{sys.argv[0]} not recognised as a Graveparser compatible program\n"\
                         f"\n sys.argv: {sys.argv}")

    return graveparser

def get_bilbypipe_args(bilby_pipe_config):
    '''
    Input:
        bilby_pipe_config - Dictionary containing settings for bilby_pipe

    Output:
        bilby_pipe_args - Known arguments to bilby_pipe
        bilby_pipe_unknown_args - Unknown arguments to bilby_pipe

    Function gets the arguments to bilby_pipe separated into known and unknown lists
    '''

    bilby_pipe_parser = create_bilbypipe_parser()
    bilby_pipe_ini = f"{bilby_pipe_config['outdir']}/bilbypipe.ini"

    argument_list = [bilby_pipe_ini]

    with open(bilby_pipe_ini, "w", encoding="utf-8") as inifile:
        for key, value in bilby_pipe_config.items():
            inifile.writelines(f"{key.replace('_','-')} = {value} \n")

    bilby_pipe_args, bilby_pipe_unknown_args = bilby_pipe_parser.parse_known_args(argument_list)

    return bilby_pipe_parser, bilby_pipe_args, bilby_pipe_unknown_args

def graveparser_lens_generation_parser():
    '''
    Output:
        graveparser - Gravelamps argument parser

    Creates an argument parser for the lens data generation utility
    '''

    prog = "gravelamps_generate_lens"
    description = "Gravelamps Lens Data Generation Utility"

    graveparser = argparse.ArgumentParser(prog=prog, description=description)

    #Position Arguments
    graveparser.add_argument('ini', metavar='INI', type=str,
                             help="Location of the INI file containing options for usage")

    #Optional Arguments Functionality
    graveparser.add_argument('-i', '--injection', action="store_true",
                             help="Will use injection lensing settings instead of analysis\
                                   settings")
    graveparser.add_argument('-l', '--local', action="store_true",
                             help="Will perform runs locally instead of generating condor submits")
    graveparser.add_argument('-s', '--submit', action="store_true",
                             help="Will directly submit condor job")

    #Optional Arguments Utility
    graveparser.add_argument('-v', '--verbose', action="store_true",
                             help="Display maximum level of information from logger")

    return graveparser

def graveparser_inference_parser():
    '''
    Output:
        graveparser - Gravelamps argument parser

    Creates an argument parser for the inference program utility
    '''

    prog = "gravelamps_inference"
    description = "Gravelamps Inference Program"

    graveparser = argparse.ArgumentParser(prog=prog, description=description)

    #Position Arguments
    graveparser.add_argument('ini', metavar='INI', type=str,
                             help="Location of the INI file containing options for usage")

    #Optional Arguments Functionality
    graveparser.add_argument('-i', '--injection', action="store_true",
                             help="Run contains injection")
    graveparser.add_argument("-l", "--local", action="store_true",
                             help="Will perform runs locally instead of generating condor submits")

    graveparser.add_argument("-s", "--submit", action="store_true",
                             help="Will submit condor jobs directly")

    #Optional Arguments Utility
    graveparser.add_argument('-v', '--verbose', action="store_true",
                             help="Display maximum level of information from logger")

    return graveparser

def graveparser_interpolator_parser():
    '''
    Output:
        graveparser - Gravelamps argument parser

    Creates an argument parser for the interpolator data generation utility for HTCondor
    '''

    prog = "gravelamps_generate_interpolator_data"
    description = "Barebones lens data generation utility. To be used by HTCondor"

    graveparser = argparse.ArgumentParser(prog=prog, description=description)

    #Position Arguments
    graveparser.add_argument('model', metavar='model', type=str,
                             help="Python path of Interpolator Generation Model")

    graveparser.add_argument('filedict', metavar='File Dict', type=str,
                              help="String containing Files for Interpolator Generation")

    graveparser.add_argument('ini', metavar='INI', type=str,
                             help="Location of the INI file containing options for usage")

    #Optional Arguments Functionality
    graveparser.add_argument('-i', '--injection', action="store_true",
                             help="Will use injection lensing settings instead of analysis\
                                   settings")
    graveparser.add_argument('-l', '--local', action="store_true",
                             help="Will perform runs locally instead of generating condor submits")

    #Optional Arguments Utility
    graveparser.add_argument('-v', '--verbose', action="store_true",
                             help="Display maximum level of information from logger")

    return graveparser
