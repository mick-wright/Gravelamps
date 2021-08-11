'''
Pipe file generators for gravelamps.inference

Functions generate the files for the pipe based inference.

Written by Mick Wright 2021
'''

import os

def lens_subfile(config, dim_freq_file, sour_pos_file, amp_fac_real_file,
                 amp_fac_imag_file, additional_lens_parameters):
    '''
    Input:
        config - INI configuration parser
        dim_freq_file - location of file containing dimensionless frequency values
        sour_pos_file - location of file containing source position values
        amp_fac_real_file - location of file containing the real part of the amplification
                            factor values
        amp_fac_imag_file - location of file containing the imaginary part of the amplification
                            factor values
        additional_lens_parameters - list containing the additional parameters for the lens model

    Function generates the condor submission file in order to generate the amplification factor
    data necessary for the analysis run
    '''

    #Get the submission subdirectory
    outdir = config.get("output_settings", "outdir")
    outdir = os.path.abspath(outdir)
    submit_directory = outdir + "/submit"

    #Get the lens model
    lens_model = config.get("lens_generation_settings", "lens_model")

    #Get the executable directory and construct the full path of the lens model function
    executable_directory = config.get(
        "lens_generation_settings", "executable_directory", fallback="default")
    if executable_directory == "default":
        executable_directory = os.path.expanduser("~") + "/bin"
    executable = executable_directory + "/" + lens_model

    #Create the submission subdirectory if it doesn't exist
    if not os.path.isdir(submit_directory):
        os.mkdir(submit_directory)

    #Create the .submit file and get the condor settings from the main INI file
    submit_file = submit_directory + "/generate_lens.sub"
    condor_settings_dictionary = config._sections["condor_settings"].copy()

    #Open the submit file and write the necessary parts
    with open(submit_file, "w") as sub:
        sub.write("universe = vanilla\n")
        sub.write("transfer_input_files = " + os.path.abspath(dim_freq_file)
                  + "/" + os.path.abspath(sour_pos_file) + "\n")
        sub.write("executable = " + executable + "\n")

        #Construct the arguments necessary for the function call
        arguments = ""

        for argument in (dim_freq_file, sour_pos_file, amp_fac_real_file, amp_fac_imag_file):
            add = os.path.abspath(argument) + " "
            arguments += add
        for argument in additional_lens_parameters:
            add = argument + " "
            arguments += add

        sub.write("arguments = " + arguments + "\n")
        sub.write("log = " + submit_directory + "/lens_generation.log\n")
        sub.write("output = " + submit_directory + "/lens_generation.out\n")
        sub.write("error = " + submit_directory + "/lens_generation.err\n")

        for key, value in condor_settings_dictionary.items():
            sub.write(key + " = " + value + "\n")

        sub.write("queue 1\n")
