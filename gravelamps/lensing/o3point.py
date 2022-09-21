'''
Isolated Point Mass Lensing Functions used within LVK O3 (full + O3a) Papers.

Written by Eungwang Seo 2020

NOTE FROM DEVELOPERS:
Included as legacy code, these functions are designed to replicate specifically those analyses,
and are superceded by Gravelamps' point mass lensing functions. Consequent to this, they will not
be actively updated. They were written to be functional for version 2.1, subsequent versions do not
guarantee that this will remain functional.

To be used, the user must specify a location for the lookuptable.h5 file that can be found:
https://git.ligo.org/eungwang.seo/bilby_pmlens/-/blob/master/src/lookuptable.h5
'''

import h5py
import numpy as np

_GCSM = 2.47701878 * 1.98892 * 10**-6
_FREQUENCY_STEP = 8.0 * np.pi * 4.93E-6 * 0.1 * 2 * 1000
_LOWER_SOURCE_STEP = 0.001
_HIGHER_SOURCE_STEP = 0.01

_higher_source_array = []
_lower_source_array = []

_lens_parameters = ["log_lens_mass", "source_position"]

def load_table(table_location):
    '''
    Input:
        table_location - string containing filepath of the lookup table of amplification factor
                         values

    Function loads in the lookup table from the specified location and generates arrays containing
    the higher and lower source position values from it.
    '''

    with h5py.File(table_location, "r") as lookup_file:
        higher_source_array_real = np.array(lookup_file.get("f11")[:])
        higher_source_array_imag = np.array(lookup_file.get("f12")[:])

        lower_source_array_real = np.array(lookup_file.get("f21")[:])
        lower_source_array_imag = np.array(lookup_file.get("f22")[:])

    higher_source_array = higher_source_array_real + 1j*higher_source_array_imag
    lower_source_array = lower_source_array_real + 1j*lower_source_array_imag

    global _higher_source_array
    global _lower_source_array

    _higher_source_array = np.transpose(higher_source_array)
    _lower_source_array = np.transpose(lower_source_array)

def magnification(source_position, mode):
    '''
    Input:
        source_position - dimensionless displacement from optical axis
        mode - flag for whether to give the positive or negative magnification

    Output:
        mag - Value of the magnification

    Function computes the mangification for the isolated point mass lens model
    '''

    numerator = 2 + source_position**2
    denominator = 2 * source_position * np.sqrt(4 + source_position**2)
    source_position_term = numerator/denominator

    if mode == 1:
        mag = 0.5 + source_position_term
    else:
        mag = 0.5 - source_position_term

    return mag

def time_delay(source_position):
    '''
    Input:
        source_position - dimensionless displacement from optical axis

    Output:
        delay - dimensionless time delay for the image

    Function computes the dimensionless time delay for the image.
    '''

    source_position_term = np.sqrt(4 + source_position**2)

    power_term = source_position * source_position_term
    log_term =\
        np.log((source_position_term + source_position)/(source_position_term - source_position))

    delay = power_term + log_term
    return delay

def frequency_to_dimensionless_frequency(frequency_array, log_lens_mass):
    '''
    Input;
        frequency_array - Array containing frequencies to be converted
        log_lens_mass - Log of the lens mass

    Ouput:
        dimensionless_frequency_array - Dimensionless equivalent array

    Function converts an array of frequnecy values to the equivalent dimensionless frequencies
    '''

    dimensionless_frequency_array = 8 * (10 ** log_lens_mass) * np.pi * frequency_array * _GCSM
    return dimensionless_frequency_array

def amplification_factor_geometric(dimensionless_frequency,
                                   source_position):
    '''
    Input:
        dimensionless_frequency - Dimensionless form of frequency being amplified
        source_position - Dimensionless displacement from the optical axis

    Output:
        amplification_factor - Value of the amplification factor

    Calculates the geometric optics amplification factor value
    '''

    mag_pos = magnification(source_position, 1)
    mag_neg = magnification(source_position, 0)
    mag_pos_term = np.sqrt(np.abs(mag_pos))
    mag_neg_term = np.sqrt(np.abs(mag_neg))

    delay = time_delay(source_position)
    delay_term = np.exp(1j * dimensionless_frequency * delay)

    factor = mag_pos_term - 1j * delay_term * mag_neg_term

    return factor

def amplification_factor_wave_lookup(dimensionless_frequency,
                                     source_position,
                                     position_index_step,
                                     data_array):
    '''
    Input:
        dimensionless_frequency - Dimensionless form of frequency being amplified
        source_position - Dimensionless displacement from the optical axis
        position_index_step - Step size for the source position
        data_array - Array of the data for the amplification factor

    Output:
        amplification_factor - Value of the amplification factor from the lookup table

    Function computes the value of the amplification factor using the lookup table.
    It does this as a correction to the geometric factor calculation
    '''

    frequency_index = int(dimensionless_frequency/_FREQUENCY_STEP)
    position_index = int((source_position - 0.1)/position_index_step)

    geometric_factor = amplification_factor_geometric(dimensionless_frequency, source_position)

    step_term = 1.0/(position_index_step * _FREQUENCY_STEP)

    position_index_base_term = (position_index + 1) * position_index_step + 0.1 - source_position
    frequency_index_base_term = (frequency_index + 1) * _FREQUENCY_STEP - dimensionless_frequency
    base_term = position_index_base_term * frequency_index_base_term\
                * data_array[position_index, frequency_index]

    frequency_index_higher_term = dimensionless_frequency - frequency_index * _FREQUENCY_STEP
    position_index_higher_term = source_position - position_index * position_index_step - 0.1

    higher_frequency_term = position_index_base_term * frequency_index_higher_term *\
                            data_array[position_index, frequency_index+1]
    higher_position_term = position_index_higher_term * frequency_index_base_term *\
                           data_array[position_index+1, frequency_index]
    higher_both_term = position_index_higher_term * frequency_index_higher_term *\
                       data_array[position_index+1, frequency_index+1]

    factor = geometric_factor + step_term\
        * (base_term + higher_frequency_term + higher_position_term + higher_both_term)

    return factor

def amplification_factor(frequency_array, source_position, log_lens_mass, **kwargs):
    '''
    Input:
        frequency_array - Array of frequencies to be amplified
        source_position - Dimensionless displacement from the optical axis
        log_lens_mass - Log of the mass of the lensing object

    Output:
        amplification_factor - Values of the amplification factor

    Function computes the values of the amplification factor from the lookup table.
    '''

    dimensionless_frequency_array =\
        frequency_to_dimensionless_frequency(frequency_array, log_lens_mass)

    amplification_factor_array = np.zeros(len(dimensionless_frequency_array))

    for idx, dimensionless_frequency in enumerate(dimensionless_frequency_array):
        if dimensionless_frequency > 246.0:
            amplification_factor_array[idx] =\
                amplification_factor_geometric(dimensionless_frequency, source_position)
        else:
            if source_position > 0.3:
                amplification_factor_array[idx] =\
                    amplification_factor_wave_lookup(dimensionless_frequency,
                                                     source_position,
                                                     _HIGHER_SOURCE_STEP,
                                                     _higher_source_array)
            else:
                amplification_factor_array[idx] =\
                    amplification_factor_wave_lookup(dimensionless_frequency,
                                                     source_position,
                                                     _LOWER_SOURCE_STEP,
                                                     _lower_source_array)

    return amplification_factor_array
