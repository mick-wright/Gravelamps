"""Gravelamps Module Handling

Following are functions that handle modules used within Gravelamps. This is largely an interface
used for handling the module that will contain the lensing functions in an agnostic fashion

Written by Mick Wright 2022

Routines
--------
get_lens_module
    Retrieves module to be used for lensing functions
check_module
    Ascertains the existence of given module
get_interpolator_model
    Retrieves module to be used for interpolator functions

"""

import importlib

from gravelamps.core.gravelog import gravelogger

def get_lens_module(config, args):
    """
    Retrieves module to be used for lensing functions

    The module is specified by the user either as a full python path or as a submodule of
    `gravelamps.lensing`. If the module is not supplied at all, it defaults to
    `gravelamps.lensing.interpolator`.

    Parameters
    ----------
    config : configparser.ConfigParser
        Object containing settings from INI file
    args : argparse.Namespace
        Object containing commandline arguments to program

    Returns
    -------
    lens_module : ModuleType
        Module containing the lensing functions
    """

    if args.injection:
        lens_type = "injection"
    else:
        lens_type = "analysis"

    lens_module = config.get(f"{lens_type}_lens_generation_settings",
                             "lensing_module", fallback=None)
    if lens_module is None:
        gravelogger.warning(
            "Lensing Module Not Specified, defaulting to gravelamps.lensing.interpolator")
        lens_module = "gravelamps.lensing.interpolator"
    elif "." not in lens_module:
        gravelogger.info("Full path not given. Assuming %s in gravelamps.lensing", lens_module)
        lens_module = f"gravelamps.lensing.{lens_module}"

    module_exists = check_module(lens_module)
    if not module_exists:
        raise ModuleNotFoundError(f"{lens_module} could not be found!")

    gravelogger.info("Using Lensing Module: %s", lens_module)

    return lens_module

def check_module(module_name):
    """
    Ascertains the existence of a given module

    Parameters
    ----------
    module_name : str
        full python path of the module to be checked

    Returns
    -------
    bool
        Flag of whether the module exists or not
    """

    module_spec = importlib.util.find_spec(module_name)
    if module_spec is None:
        return False
    return True

def get_interpolator_model(config, args):
    """
    Retrieves module to be used for interpolator functions

    The model for the interpolator should be specified in the INI as either a complete python path
    or one of the submodules of `gravelamps.lensing`. Will return None if no module is specified

    Parameters
    ----------
    config : configparser.ConfigParser
        Object containing settings from the INI file
    args : argparse.Namespace
        Object containing commandline arguments to the program

    Returns
    -------
    interpolator_model : str
        Full python path to the interpolator model module

    Raises
    ------
    ModuleNotFoundError
        If the full python path as determined by the function cannot be loaded as a module
    """

    if args.injection:
        lens_type = "injection"
    else:
        lens_type = "analysis"

    interpolator_model =\
        config.get(f"{lens_type}_lens_generation_settings", "interpolator_model", fallback=None)

    if interpolator_model is None:
        gravelogger.warning(
            "No Interpolator Model specified, cannot generate new interpolator data!")
        return None

    if "." in interpolator_model:
        if not check_module(interpolator_model):
            raise ModuleNotFoundError(f"Interpolator model {interpolator_model} not found!")
    else:
        gravelogger.info("Assuming %s in gravelamps.lensing", interpolator_model)
        interpolator_model = f"gravelamps.lensing.{interpolator_model}"
        if not check_module(interpolator_model):
            raise ModuleNotFoundError(f"Interpolator model {interpolator_model} not found!")

    gravelogger.info("Using Interpolator Model: %s", interpolator_model)
    return interpolator_model
