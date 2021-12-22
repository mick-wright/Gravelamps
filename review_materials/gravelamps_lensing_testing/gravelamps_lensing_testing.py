'''
gravelamps.lensing Testing

The script below uses a minimal example INI to test the functionality of the
classes and functions contained within gravelamps.lensing

Written by Mick Wright 2021
'''

import os
from configparser import ConfigParser

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import bilby
from astropy.constants import M_sun

import gravelamps.lensing
import gravelamps.inference

matplotlib.rcParams.update({"text.usetex":True, "font.family":"serif"})

#Read in INI
config = ConfigParser()
config.read("lensing_testing.ini")

#Set up output folders
outdir = config.get("output_settings", "outdir")

if not os.path.isdir(outdir):
    os.mkdir(outdir)
    os.mkdir(outdir+"/data")

#Read in interpolator files
dim_freq_file = config.get("lens_generation_settings", "dimensionless_frequency_file")
sour_pos_file = config.get("lens_generation_settings", "source_position_file")
amp_fac_real_file = config.get("lens_generation_settings", "amplification_factor_real_file")
amp_fac_imag_file = config.get("lens_generation_settings", "amplification_factor_imag_file")

#Get the Waveform Generator and Frequency Domain Source Model
waveform_generator_class = config.get("analysis_settings", "lensed_waveform_generator_class")
frequency_domain_source_model = config.get(
        "analysis_settings", "lensed_frequency_domain_source_model")

waveform_generator_class, frequency_domain_source_model = (
        gravelamps.inference.helpers.wfgen_fd_source(
            waveform_generator_class, frequency_domain_source_model))

#Set up the Waveform Arguments dictionary
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

#Read in Duration and Sampling Frequency
duration = config.getfloat("analysis_settings", "duration")
sampling_frequency = config.getfloat("analysis_settings", "sampling_frequency")

#Set up the Waveform Generator
waveform_generator = waveform_generator_class(
        duration=duration, sampling_frequency=sampling_frequency,
        frequency_domain_source_model=frequency_domain_source_model,
        parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
        waveform_arguments=waveform_arguments)

#Set up the Interferometers and Inject Signal
interferometer_list = config.get(
        "analysis_settings", "interferometers").replace(" ","").split(",")

interferometers = bilby.gw.detector.InterferometerList(interferometer_list)

injection_parameters = config._sections["injection_parameters"].copy()
injection_parameters.update(
        (key,float(value)) for key, value in injection_parameters.items())

trigger_time = config.getfloat("analysis_settings", "trigger_time")
gps_start_time = trigger_time + 2 - duration

injection_parameters["geocent_time"] = trigger_time

strains = waveform_generator.frequency_domain_strain(injection_parameters)

interferometers.set_strain_data_from_power_spectral_densities(
        sampling_frequency=sampling_frequency, duration=duration, start_time=gps_start_time)
interferometers.inject_signal(
        waveform_generator=waveform_generator, parameters=injection_parameters)

interferometers.plot_data(label="waveform_testing")

plt.figure(figsize=(16,10))
plt.plot(interferometers[0].frequency_array, np.abs(strains["plus"]), label="$|h_{+}|$")
plt.plot(interferometers[0].frequency_array, np.abs(strains["cross"]), label=r"$|h_{\times}|$")
plt.title("Waveform Generator Strains for 36/29 Object as Lensed by 50 Solar Mass object")
plt.legend()
plt.savefig("waveform-generator-strains.pdf", bbox_inches="tight")

#Simple testing of the dimensionless frequency and natural mass converter functions
lens_mass = 50
frequency_test_value = 600

lens_mass_kg = 50 * M_sun.value

nat_mass = gravelamps.lensing.natural_mass(lens_mass)
nat_mass_kg = gravelamps.lensing.natural_mass(lens_mass_kg, mode="kg")

print(f"A Mass of {lens_mass} Solar Masses is in natural units: {nat_mass}")
print(f"A mass of {lens_mass_kg} kg (same as above) is in natural units: {nat_mass_kg}")

dim_freq = gravelamps.lensing.dimensionless_frequency(frequency_test_value, nat_mass)

print(f"The dimensionless frequency equivalent of {frequency_test_value} Hz lensed by that mass is {dim_freq}")
