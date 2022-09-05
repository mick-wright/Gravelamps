'''
Gravelamps Utility Functions

Functions within are utility functions for the operation of Gravelamps. They are unliekely to be
necessary for any use by a user.

Written by Mick Wright 2022
'''

import importlib

from gravelamps.core.gravelog import gravelogger

def get_lens_module(config, args):
    '''
    Input:
        config - INI configuration parser
        args - Commandline arguments given to program

    Output:
        lens_module - Name of the Module containing the lensing functions

    Function checks the user configuraiton for the lensing module, defaulting to
    gravelamps.lensing.interpolator and verifies that the module can be loaded.
    '''

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
    '''
    Input:
        module_name - Name of model to check existence of

    Function checks existence of module given, returning a boolean based on that
    '''

    module_spec = importlib.util.find_spec(module_name)
    if module_spec is None:
        return False
    return True

def get_interpolator_model(config, args):
    '''
    Input:
        config - INI configuration parser
        args - Commandline arguments given to program

    Output:
        interpolator_model - Model specified by user defaulting to None

    Function will retrieve the interpolator model from the INI and will check that either it is a
    module in it's own right or is a submodule of gravelamps.lensing. If neither will return None
    '''

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
