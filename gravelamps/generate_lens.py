'''
Gravelamps Lens Data Generation

Functions within control generating lensing data.
The main function will generate lensing data for the interpolaor methods.

Written by Mick Wright 2022
'''

import os

from gravelamps.core.gravelog import gravelogger, setup_logger
from gravelamps.core.graveparser import create_graveparser
from gravelamps.core.file_handling import (get_config,
                                         retrieve_interpolator_files,
                                         grid_file_handler,
                                         get_output_directories,
                                         data_file_handler)
from gravelamps.core.module_handling import get_lens_module, get_interpolator_model
from gravelamps.lensing.generic import (generate_interpolator_data,
                                        generate_interpolator_condor,
                                        get_condor_config)

def main(_config=None, _args=None, _injection=False):
    '''
    Optional Input:
        _config - INI configuration parser
        _args - Commandline arguments from program
        _injection - Flag to set args.injection

    Ouptut:
        interpolator_files - Dictionary of files containing interpolator data

    Function will use input parameters from the user defined INI to generate lensing data
    for use with the interpolator methods
    '''

    if _config is None:
        graveparser = create_graveparser()
        args = graveparser.parse_args()

        config = get_config(args)

        output_directories = get_output_directories(config)
        logging_level = config.get("output_settings", "logging_level", fallback="INFO")
        setup_logger(output_directories["outdir"], logging_level, args=args)
        gravelogger.info("Gravelogger setup complete, log will be output to %s/grave.log",
                         output_directories["outdir"])

        if config.getboolean("run_settings", "injection", fallback=False):
            args.injection = True
        if config.getboolean("run_settings", "local", fallback=False):
            args.local = True
        if config.getboolean("run_settings", "submit", fallback=False):
            args.submit = True

    else:
        config = _config
        args = _args
        args.injection = _injection

        output_directories = get_output_directories(config)

    if args.injection:
        gravelogger.info("Running in Injection mode")

    lens_module = get_lens_module(config, args)
    if lens_module == "gravelamps.lensing.interpolator":
        interpolator_model = get_interpolator_model(config, args)
    else:
        interpolator_model = lens_module

    interpolator_files = retrieve_interpolator_files(config, args)
    interpolator_files = grid_file_handler(config, args,
                                           output_directories["data"],
                                           interpolator_files)
    interpolator_files, files_complete = data_file_handler(args,
                                                           output_directories["data"],
                                                           interpolator_files)

    if files_complete < 2:
        if interpolator_model is None:
            raise IOError(("Amplification factor data incomplete with no defined interpolator"\
                           " model. Cannot create new data to complete interpolator!"))
        gravelogger.info("Amplification factor data incomplete, will use %s to generate",
                         interpolator_model)

        if args.local:
            generate_interpolator_data(config, args, interpolator_model, interpolator_files)
        else:
            condor_settings = get_condor_config(config, args,
                                                output_directories,
                                                interpolator_model,
                                                interpolator_files)
            generate_interpolator_condor(args, output_directories, condor_settings)

    if _config is not None:
        return interpolator_files
    return

if __name__ == "__main__":
    main()
