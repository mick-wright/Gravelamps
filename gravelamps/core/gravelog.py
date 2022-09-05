'''
Gravelamps Logging

Functions within control the use of the logger used within Gravelamps utilising the python builtin
logging module

Written by Mick Wright 2022
'''

import logging
import os

def setup_logger(outdir=None, log_level="INFO", args=None):
    '''
    Input:
        outdir - top level directory containing run output
        log_level - string containing either the text or numerical logging level
        args - command line arguments given to program

    Function generates the logger to be used throughout Gravelamps
    '''

    #Logging level is overridden with highest possible if verbose output is requested
    if args is not None:
        if args.verbose:
            log_level = "DEBUG"

    #For a string log level, get the level from the logging module
    if isinstance(log_level, str):
        try:
            log_level = getattr(logging, log_level.upper())
        except AttributeError as exc:
            raise ValueError(f"{log_level.upper()} is not a valid logging level!") from exc
    else:
        log_level = int(log_level)

    logger = logging.getLogger("Gravelamps")
    logger.handlers.clear()
    logger.propagate = False
    logger.setLevel(log_level)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    logger.addHandler(console_handler)

    if outdir is not None:
        if not os.path.isdir(outdir):
            os.mkdir(outdir)

        log_file = f"{outdir}/grave.log"
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(log_level)

        formatter = logging.Formatter("%(asctime)s %(name)s %(levelname)s: %(message)s")
        formatter.datefmt = "%H:%M:%S %d/%m/%Y"
        console_handler.setFormatter(formatter)
        file_handler.setFormatter(formatter)

        logger.addHandler(file_handler)

    logger.debug("Logger activated")

setup_logger()
gravelogger = logging.getLogger("Gravelamps")
