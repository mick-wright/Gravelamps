'''
Testing gravelamps.lensing.utils

The script below uses simple testing cases to test the output of the
gravelamps.lensing.utils file

Written by Mick Wright 2021
'''

from configparser import ConfigParser
import os
import numpy as np
import matplotlib.pyplot as plt
import gravelamps.lensing

#Read in the simple INI file utils_testing.ini
config = ConfigParser()
config.read("utils_testing.ini")

#Check if the data directory exists, if it doesn't make it
outdir = config.get("output_settings", "outdir")
data_subdirectory = outdir + "/data"

if not os.path.isdir(data_subdirectory):
    os.mkdir(data_subdirectory)

#Test generation of the dimensionless frequency file
w_file = gravelamps.lensing.utils.generate_dimensionless_frequency_file(config)
print(f"Dimensionless Frequency file located at {w_file}")

#Test generation of the source position file
y_file = gravelamps.lensing.utils.generate_source_position_file(config)
print(f"Source Position file located at {y_file}")

#Test the generation of the interpolator

#Generate data covering dimensionless frequency and source position
dim_freq_array = np.linspace(0, 2*np.pi, 120)
sour_pos_array = np.linspace(0, 2*np.pi, 100)

#Generating amplification factor data - using dummy sin functions 
amp_fac_real_array = np.zeros((len(dim_freq_array), len(sour_pos_array)))
amp_fac_imag_array = np.zeros(amp_fac_real_array.shape)

for i in range(len(sour_pos_array)):
    amp_fac_real_array[:,i] = np.sin(dim_freq_array)
for j in range(len(dim_freq_array)):
    amp_fac_imag_array[j,:] = np.sin(sour_pos_array)

#Create files containing test data
np.savetxt("test_w.dat", dim_freq_array)
np.savetxt("test_y.dat", sour_pos_array)
np.savetxt("test_freal.dat", amp_fac_real_array)
np.savetxt("test_fimag.dat", amp_fac_imag_array)

#Generate interpolator
interpolator = gravelamps.lensing.utils.generate_interpolator(
    "test_w.dat", "test_y.dat", "test_freal.dat", "test_fimag.dat")

#Shift arrays to test interpolation 
shifted_dim_freq = dim_freq_array + dim_freq_array[1]/2
shifted_sour_pos = sour_pos_array + sour_pos_array[1]/2

#Plots show the dimensionless frequency interpolation 
#Left: Calculated values vs Interpolated 
#Right: Residual between the two
#Top: Using the exact data points used to generate interpolator
#Bottom: Shifting the data to use non of the points used to generate interpolator
plt.subplot(2,2,1)
plt.plot(dim_freq_array, np.sin(dim_freq_array), label="Calculated")
plt.plot(dim_freq_array, np.real(interpolator(dim_freq_array, 0)),
         linestyle='--', color='k', label="Interpolated")
plt.legend()
plt.subplot(2,2,2)
plt.plot(dim_freq_array, np.real(interpolator(dim_freq_array, 0)) - np.sin(dim_freq_array))
plt.subplot(2,2,3)
plt.plot(shifted_dim_freq, np.sin(shifted_dim_freq), label="Calculated")
plt.plot(shifted_dim_freq, np.real(interpolator(shifted_dim_freq, 0)),
         linestyle="--", color='k', label="Interpolated")
plt.legend()
plt.subplot(2,2,4)
plt.plot(shifted_dim_freq, np.real(interpolator(shifted_dim_freq, 0)) - np.sin(shifted_dim_freq))
plt.show()

#Same as above but for source position interpolation
plt.subplot(2,2,1)
plt.plot(sour_pos_array, np.sin(sour_pos_array), label="Caluclated")
plt.plot(sour_pos_array, np.imag(interpolator(0, sour_pos_array)),
         linestyle="--", color='k', label="Interpolated")
plt.legend()
plt.subplot(2,2,2)
plt.plot(sour_pos_array, np.imag(interpolator(0, sour_pos_array)) - np.sin(sour_pos_array))
plt.subplot(2,2,3)
plt.plot(shifted_sour_pos, np.sin(shifted_sour_pos),label="Calculated")
plt.plot(shifted_sour_pos, np.imag(interpolator(0, shifted_sour_pos)),
         linestyle="--", color='k', label="Interpolated")
plt.legend()
plt.subplot(2,2,4)
plt.plot(shifted_sour_pos, np.imag(interpolator(0, shifted_sour_pos)) - np.sin(shifted_sour_pos))
plt.show()
