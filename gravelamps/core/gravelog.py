"""Gravelamps Logging

Functions within control the use and generation of the gravelogger---implementation of the logging
module in use within gravelamps.

Routines
--------
setup_logger
    Generates/modifies the gravelogger used throughout Gravelamps functions

Notes
-----
Importing the module generates a blank, default gravelogger which is then modified whenever
`setup_logger` is called. By default, this will not add text logging to avoid spam of grave.log
output files. 

Written by Mick Wright 2022
"""

import logging
import os

def setup_logger(outdir=None, log_level="INFO", args=None):
    """
    Generates/modifies the gravelogger used throughout Gravelamps functions

    By default, the logger will simply output to the interpreter. If an output directory is
    specified, the logger will also save output to grave.log within that directory. Logging
    levels may be specified.

    Parameters
    ----------
    outdir : str, optional
        Path in which to generate a grave.log output saving the messages logged
    log_level: str, optional, default='INFO'
        Level of logging, can be `INFO`, `WARNING`, `ERROR`, `DEBUG`
    args : argparse.Namespace, optional
        Object containing commandline arguments to program
    """

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
