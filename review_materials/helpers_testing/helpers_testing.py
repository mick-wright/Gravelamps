'''
gravelamps.inference.helpers Testing

The script below uses minimal example INIs to test the functionality of those
functions inside gravelamps.inference.helpers

Written by Mick Wright 2021
'''

import os
import shutil
from configparser import ConfigParser

import gravelamps.inference.helpers

#Read in the four different configurations
config_files_no_copy = ConfigParser()
config_files_w_copy = ConfigParser()
config_no_files = ConfigParser()
config_non_standard = ConfigParser()

config_files_no_copy.read("files_no_copy.ini")
config_files_w_copy.read("files_w_copy.ini")
config_no_files.read("no_files.ini")
config_non_standard.read("non_standard_model.ini")

if os.path.isdir(config_no_files.get("output_settings", "outdir")):
    shutil.rmtree(config_no_files.get("output_settings", "outdir"))

for config in [config_files_no_copy, config_files_w_copy, config_no_files, config_non_standard]:
    if not os.path.isdir(config.get("output_settings", "outdir")):
        os.mkdir(config.get("output_settings", "outdir"))
        os.mkdir(config.get("output_settings", "outdir")+"/data")
        os.mkdir(config.get("output_settings", "outdir")+"/submit")

#wyhandler testing for each of the file cases
w_fnc, y_fnc = gravelamps.inference.helpers.wy_handler(config_files_no_copy)
w_fwc, y_fwc = gravelamps.inference.helpers.wy_handler(config_files_w_copy)
w_nf, y_nf = gravelamps.inference.helpers.wy_handler(config_no_files)
w_ns, y_ns = gravelamps.inference.helpers.wy_handler(config_non_standard)

print(f"For files_no_copy, w file is located at {w_fnc} and y file at {y_fnc}")
print(f"For files_w_copy, w file is located at {w_fwc} and y file at {y_fwc}")
print(f"For no_files, w_file is located at {w_nf} and y file at {y_nf}")

#Additional parameter testing
pointlens_parameters = gravelamps.inference.helpers.get_additional_parameters(config_no_files)
sis_parameters = gravelamps.inference.helpers.get_additional_parameters(config_files_no_copy)
nfw_parameters = gravelamps.inference.helpers.get_additional_parameters(config_files_w_copy)
ns_parameters = gravelamps.inference.helpers.get_additional_parameters(config_non_standard)

print(f"Additional parameters for point mass: {pointlens_parameters}")
print(f"Additional parameters for sis: {sis_parameters}")
print(f"Additional parameters for nfw: {nfw_parameters}")
print(f"Additional parameters for non-standard parameters: {ns_parameters}")

#amp_fac_handler testing
fr_fnc, fi_fnc = gravelamps.inference.helpers.amp_fac_handler(
    config_files_no_copy, w_fnc, y_fnc)
fr_fwc, fi_fwc = gravelamps.inference.helpers.amp_fac_handler(
    config_files_w_copy, w_fwc, y_fwc)
fr_nf, fi_nf = gravelamps.inference.helpers.amp_fac_handler(
    config_no_files, w_nf, y_nf)
fr_nfp, fi_nfp = gravelamps.inference.helpers.amp_fac_handler(
    config_non_standard, w_nf, y_nf, mode="pipe")

print(f"For files_no_copy. real file located at {fr_fnc}, imag at {fi_fnc}")
print(f"For files_w_copy, real file located at {fr_fwc}, imag at {fi_fwc}")
print(f"For no_files, read file located at {fr_fwc}, imag at {fi_fwc}")

#wfgen_fd_source testing
wfgen_split_example = "bilby.gw.waveform_generator.WaveformGenerator"
wfgen_default_example = "default"
wfgen_single_name_example = "WaveformGenerator"
wfgen_nonsense_example = "nonsense.example"

fds_split_example = "bilby.gw.source.lal_binary_black_hole"
fds_single_example = "lal_binary_black_hole"
fds_nonsense_example = "nonsense.example"

wfgen_split, fds_split = gravelamps.inference.helpers.wfgen_fd_source(
    wfgen_split_example, fds_split_example)
print(wfgen_split)
print(fds_split)
wfgen_default, fds_default = gravelamps.inference.helpers.wfgen_fd_source(
    wfgen_default_example, fds_split_example)
print(wfgen_default)
wfgen_single_name, fds_single = gravelamps.inference.helpers.wfgen_fd_source(
    wfgen_single_name_example, fds_single_example)
print(wfgen_single_name)
print(fds_single)

try:
    wfgen_nonsense, fds_nonsense = gravelamps.inference.helpers.wfgen_fd_source(
        wfgen_nonsense_example, fds_nonsense_example)
except:
    print("Nonsense Waveform Generator Failed")

try:
    wfgen_nonsense, fds_nonsense = gravelamps.inference.helpers.wfgen_fd_source(
        wfgen_single_name_example, fds_nonsense_example)
except:
    print("Nonsense Frequency Domain Source Model Failed")
