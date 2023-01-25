"""Gravelamps Program Argument Parsing

Following are functions handling the various argument parsers used by programs within Gravelamps
using the python builtin argparse module.

Written by Mick Wright 2022

Routines
--------
create_bilbypipe_parser
    Create a parser for bilby_pipe function usage
create_graveparser
    Metafunction that creates the specific gravelamps parser based on the program name
get_bilbypipe_args
    Parses arguments into known and unknown for bilby_pipe
graveparser_lens_generation_parser
    Creates a parser for lens generation program usage
graveparser_inference_parser
    Creates a parser for inference program usage
graveparser_interpolator_parser
    Creates a parser for the minimal interpolator data program usage

"""

import argparse
import json
import os
import sys

from bilby_pipe.parser import create_parser

def create_bilbypipe_parser():
    """
    Creates a parser for bilby_pipe function usage

    Returns
    -------
    bilbypipe_parser : argparse.ArgumentParser
        parser for the bilby_pipe program

    See Also
    --------
    bilby_pipe documentation: for a complete list of arguments that can be parsed by this parser
    """

    bilbypipe_parser = create_parser()

    return bilbypipe_parser

def create_graveparser():
    """
    Metafunction that creates the specific gravelamps parser based on the program name

    Returns
    -------
    graveparser : argparse.ArgumentParser
        parser for the specific gravelamps program being run
    """

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
    """
    Parses a configuration file for bilby_pipe into known and unknown arguments

    Parameters
    ----------
    bilby_pipe_config : dict
        Contains arguments to bilby_pipe to be stored in configuration file

    Returns
    -------
    bilby_pipe_parser : argparse.ArgumentParser
        bilby_pipe configuration parser
    bilby_pipe_args : argparse.Namespace
        Object containing known settings and values for bilby_pipe
    bilby_pipe_unknown_args : argparse.Namespace
        Object containing any unknown settings and values for bilby_pipe
    """

    bilby_pipe_parser = create_bilbypipe_parser()
    bilby_pipe_ini = f"{bilby_pipe_config['outdir']}/bilbypipe.ini"

    argument_list = [bilby_pipe_ini]

    with open(bilby_pipe_ini, "w", encoding="utf-8") as inifile:
        for key, value in bilby_pipe_config.items():
            inifile.writelines(f"{key.replace('_','-')} = {value} \n")

    bilby_pipe_args, bilby_pipe_unknown_args = bilby_pipe_parser.parse_known_args(argument_list)

    return bilby_pipe_parser, bilby_pipe_args, bilby_pipe_unknown_args

def graveparser_lens_generation_parser():
    """
    Create a parser for lens generation program usage

    Returns
    -------
    graveparser : argparse.ArgumentParser
        `gravelamps_generate_lens` program argument parser
    """

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
    """
    Create a parser for inference program usage

    Returns
    -------
    graveparser : argparse.ArgumentParser
        `gravelamps_inference` program argument parser

    """

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
    """
    Create a parser for the minimal lens interpolator program usage

    Returns
    -------
    graveparser : argparse.ArgumentParser
        Minimal lens interpolator (`gravelamps_interpolator_data`) program argument parser
    """

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
