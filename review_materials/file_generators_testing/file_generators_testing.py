'''
gravelamps.inference.file_generators Testing

The script below uses minimal example INIs to test the functionalities of those
functions inside gravelamps.inference.file_generators

Written by Mick Wright 2021
'''

import os
from configparser import ConfigParser

import gravelamps.inference

#Read in the INI files
config_lensed = ConfigParser()
config_unlensed = ConfigParser()

config_lensed.read("file_generators_lensed_testing.ini")
config_unlensed.read("file_generators_unlensed_testing.ini")

#Make all of the necessary output directories
for config in [config_lensed, config_unlensed]:
    if not os.path.isdir(config.get("output_settings", "outdir")):
        os.mkdir(config.get("output_settings", "outdir"))
        os.mkdir(config.get("output_settings", "outdir")+"/data")
        os.mkdir(config.get("output_settings", "outdir")+"/submit")

#Load in the data files and gather the additional parameters
dim_freq_file = config.get("lens_generation_settings", "dimensionless_frequency_file")
sour_pos_file = config.get("lens_generation_settings", "source_position_file")
amp_fac_real_file = config.get("lens_generation_settings", "amplification_factor_real_file")
amp_fac_imag_file = config.get("lens_generation_settings", "amplification_factor_imag_file")
additional_lens_parameters = gravelamps.inference.helpers.get_additional_parameters(config_lensed)

#Generate lens subfile
gravelamps.inference.file_generators.lens_subfile(
    config_lensed, dim_freq_file, sour_pos_file, amp_fac_real_file,
    amp_fac_imag_file, additional_lens_parameters)

#Load in Injection Parameters
injection_parameters = config_lensed._sections["injection_parameters"].copy()
injection_parameters.update(
    (key, float(value)) for key, value in injection_parameters.items())

#Generation Injection FIle
inject_file = gravelamps.inference.file_generators.injection_file(
    config, injection_parameters)

print(f"Injection file located at: {inject_file}")

#Generate bilby pipe INI file
waveform_arguments = {}

waveform_approximant = config.get("analysis_settings", "waveform_approximant")
minimum_frequency = config.getfloat("analysis_settings", "minimum_frequency")
maximum_frequency = config.getfloat("analysis_settings", "maximum_frequency")
reference_frequency = config.getfloat("analysis_settings", "reference_frequency")

waveform_arguments["waveform_approximant"] = waveform_approximant
waveform_arguments["minimum_frequency"] = minimum_frequency
waveform_arguments["maximum_frequency"] = maximum_frequency
waveform_arguments["reference_frequency"] = reference_frequency
waveform_arguments["dim_freq_file"] = dim_freq_file
waveform_arguments["sour_pos_file"] = sour_pos_file
waveform_arguments["amp_fac_real_file"] = amp_fac_real_file
waveform_arguments["amp_fac_imag_file"] = amp_fac_imag_file

injection_waveform_arguments = waveform_arguments.copy()

#Generate lensed and unlensed bilby pipe INI files
lensed_ini = gravelamps.inference.file_generators.bilby_pipe_ini(
    config=config_lensed, inject_file=inject_file,
    injection_waveform_arguments=injection_waveform_arguments,
    waveform_arguments=waveform_arguments, mode="lensed")

unlensed_ini = gravelamps.inference.file_generators.bilby_pipe_ini(
    config=config_unlensed, inject_file=inject_file,
    injection_waveform_arguments=injection_waveform_arguments,
    waveform_arguments=waveform_arguments, mode="unlensed")

print(f"Lensed only INI located at: {lensed_ini}")
print(f"Unlensed + Lensed INI located at: {unlensed_ini}")

#Finally generate an overarching DAG file - using unlensed case to show larger file
unlensed_dag = gravelamps.inference.file_generators.overarching_dag(config_unlensed)

print(f"Unlensed + Lensed DAG file located at: {unlensed_dag}")
